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
Include file that defines (possibly multiple) systems of equations.
</div>

<div class ="fragment">
<pre>
        #include "equation_systems.h"
</pre>
</div>
<div class = "comment">
Include file that defines a simple steady system
</div>

<div class ="fragment">
<pre>
        #include "transient_system.h"
        
        int main (int argc, char** argv)
        {
          libMesh::init (argc, argv);
</pre>
</div>
<div class = "comment">
This set of braces are used to force object scope.  This
way we can guarantee all our objects are destroyed before calling
libMesh::close() at the end of the program.
</div>

<div class ="fragment">
<pre>
          {
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
Create a 2D mesh.
</div>

<div class ="fragment">
<pre>
            Mesh mesh (2);
            
</pre>
</div>
<div class = "comment">
Use the internal mesh generator to create a uniform
grid on the unit square.  By default a mesh of QUAD4
elements will be created.  We instruct the mesh generator
to build a mesh of 5x5 elements.
</div>

<div class ="fragment">
<pre>
            mesh.build_cube (5, 5);
            
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
            equation_systems.set_flag ("test");
              
</pre>
</div>
<div class = "comment">
Set a simulation-specific parameter visible for all systems.
This helps in inter-system-communication.
</div>

<div class ="fragment">
<pre>
            equation_systems.set_parameter ("dummy") = 42.;
              
</pre>
</div>
<div class = "comment">
Set another simulation-specific parameter 
</div>

<div class ="fragment">
<pre>
            equation_systems.set_parameter ("nobody") = 0.;
            
</pre>
</div>
<div class = "comment">
Now we declare the system and its variables.
We begin by adding a "SteadyStytem" to the
EquationSystems object, and we give it the name
"Simple System".
</div>

<div class ="fragment">
<pre>
            equation_systems.add_system&lt;TransientImplicitSystem&gt; ("Simple System");
              
</pre>
</div>
<div class = "comment">
Adds the variable "u" to "Simple System".  "u"
will be approximated using first-order approximation.
</div>

<div class ="fragment">
<pre>
            equation_systems("Simple System").add_variable("u", FIRST);
              
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
            if (argc == 2)
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
All our objects will be destroyed at this closing brace.
That way we can safely call PetscFinalize() and be sure we
don't have any Petsc-dependent objects lurking around!
</div>

<div class ="fragment">
<pre>
          }
        
</pre>
</div>
<div class = "comment">
Call libMesh::close(), which in turn finalizes PETSc and/or
MPI.
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
  #include <FONT COLOR="#BC8F8F"><B>&quot;equation_systems.h&quot;</FONT></B>
  #include <FONT COLOR="#BC8F8F"><B>&quot;transient_system.h&quot;</FONT></B>
  
  <FONT COLOR="#228B22"><B>int</FONT></B> main (<FONT COLOR="#228B22"><B>int</FONT></B> argc, <FONT COLOR="#228B22"><B>char</FONT></B>** argv)
  {
    libMesh::init (argc, argv);
    {
      std::cout &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot;Running &quot;</FONT></B> &lt;&lt; argv[0];
      <B><FONT COLOR="#A020F0">for</FONT></B> (<FONT COLOR="#228B22"><B>int</FONT></B> i=1; i&lt;argc; i++)
        std::cout &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot; &quot;</FONT></B> &lt;&lt; argv[i];
      std::cout &lt;&lt; std::endl &lt;&lt; std::endl;
      
      Mesh mesh (2);
      
      mesh.build_cube (5, 5);
      
      mesh.print_info();
  
      EquationSystems equation_systems (mesh);
      
      equation_systems.set_flag (<FONT COLOR="#BC8F8F"><B>&quot;test&quot;</FONT></B>);
        
      equation_systems.set_parameter (<FONT COLOR="#BC8F8F"><B>&quot;dummy&quot;</FONT></B>) = 42.;
        
      equation_systems.set_parameter (<FONT COLOR="#BC8F8F"><B>&quot;nobody&quot;</FONT></B>) = 0.;
      
      equation_systems.add_system&lt;TransientImplicitSystem&gt; (<FONT COLOR="#BC8F8F"><B>&quot;Simple System&quot;</FONT></B>);
        
      equation_systems(<FONT COLOR="#BC8F8F"><B>&quot;Simple System&quot;</FONT></B>).add_variable(<FONT COLOR="#BC8F8F"><B>&quot;u&quot;</FONT></B>, FIRST);
        
      equation_systems.init();
        
      equation_systems.print_info();
  
      <B><FONT COLOR="#A020F0">if</FONT></B> (argc == 2)
        <B><FONT COLOR="#A020F0">if</FONT></B> (argv[1][0] != <FONT COLOR="#BC8F8F"><B>'-'</FONT></B>)
  	{
  	  std::cout &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot;&lt;&lt;&lt; Writing system to file &quot;</FONT></B> &lt;&lt; argv[1]
  		    &lt;&lt; std::endl;
  	  
  	  equation_systems.write (argv[1], libMeshEnums::WRITE);
  	  
  	  equation_systems.clear ();
  
  	  std::cout &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot;&gt;&gt;&gt; Reading system from file &quot;</FONT></B> &lt;&lt; argv[1]
  		    &lt;&lt; std::endl &lt;&lt; std::endl;
  	  
  	  equation_systems.read (argv[1], libMeshEnums::READ);
  
  	  equation_systems.print_info();
  	}
      
    }
  
    <B><FONT COLOR="#A020F0">return</FONT></B> libMesh::close();
  }
</pre> 
<a name="output"></a> 
<br><br><br> <h1> The console output of the program: </h1> 
<pre>
Linking ex2...
/home/peterson/code/libmesh/contrib/tecplot/lib/i686-pc-linux-gnu/tecio.a(tecxxx.o)(.text+0x1a7): In function `tecini':
: the use of `mktemp' is dangerous, better use `mkstemp'
***************************************************************
* Running Example  ./ex2
***************************************************************
 
Running ./ex2

 Mesh Information:
  mesh_dimension()=2
  spatial_dimension()=3
  n_nodes()=36
  n_elem()=25
   n_local_elem()=25
   n_active_elem()=25
  n_subdomains()=1
  n_processors()=1
  processor_id()=0

 EquationSystems
  n_systems()=1
   System "Simple System"
    Type "TransientImplicit"
    Variables="u" 
    Finite Element Types="0" 
    Approximation Orders="1" 
    n_dofs()=36
    n_local_dofs()=36
    n_constrained_dofs()=0
    n_vectors()=1
  n_flags()=1
   Flags:
    "test"
  n_parameters()=4
   Parameters:
    "dummy"=42
    "linear solver maximum iterations"=5000
    "linear solver tolerance"=1e-12
    "nobody"=0


 ---------------------------------------------------------------------------- 
| Reference count information                                                |
 ---------------------------------------------------------------------------- 
| 12SparseMatrixIdE reference count information:
| Creations:    1
| Destructions: 1
| 13NumericVectorIdE reference count information:
| Creations:    5
| Destructions: 5
| 21LinearSolverInterfaceIdE reference count information:
| Creations:    1
| Destructions: 1
| 4Elem reference count information:
| Creations:    125
| Destructions: 125
| 4Node reference count information:
| Creations:    36
| Destructions: 36
| 6DofMap reference count information:
| Creations:    1
| Destructions: 1
| 6System reference count information:
| Creations:    1
| Destructions: 1
 ---------------------------------------------------------------------------- 
 
Running ./ex2 eqn_sys.dat

 Mesh Information:
  mesh_dimension()=2
  spatial_dimension()=3
  n_nodes()=36
  n_elem()=25
   n_local_elem()=25
   n_active_elem()=25
  n_subdomains()=1
  n_processors()=1
  processor_id()=0

 EquationSystems
  n_systems()=1
   System "Simple System"
    Type "TransientImplicit"
    Variables="u" 
    Finite Element Types="0" 
    Approximation Orders="1" 
    n_dofs()=36
    n_local_dofs()=36
    n_constrained_dofs()=0
    n_vectors()=1
  n_flags()=1
   Flags:
    "test"
  n_parameters()=4
   Parameters:
    "dummy"=42
    "linear solver maximum iterations"=5000
    "linear solver tolerance"=1e-12
    "nobody"=0

<<< Writing system to file eqn_sys.dat
>>> Reading system from file eqn_sys.dat

 EquationSystems
  n_systems()=1
   System "Simple System"
    Type "TransientImplicit"
    Variables="u" 
    Finite Element Types="0" 
    Approximation Orders="1" 
    n_dofs()=36
    n_local_dofs()=36
    n_constrained_dofs()=0
    n_vectors()=1
  n_flags()=1
   Flags:
    "test"
  n_parameters()=4
   Parameters:
    "dummy"=42
    "linear solver maximum iterations"=5000
    "linear solver tolerance"=1e-12
    "nobody"=0


 ---------------------------------------------------------------------------- 
| Reference count information                                                |
 ---------------------------------------------------------------------------- 
| 12SparseMatrixIdE reference count information:
| Creations:    2
| Destructions: 2
| 13NumericVectorIdE reference count information:
| Creations:    10
| Destructions: 10
| 21LinearSolverInterfaceIdE reference count information:
| Creations:    2
| Destructions: 2
| 4Elem reference count information:
| Creations:    125
| Destructions: 125
| 4Node reference count information:
| Creations:    36
| Destructions: 36
| 6DofMap reference count information:
| Creations:    2
| Destructions: 2
| 6System reference count information:
| Creations:    2
| Destructions: 2
 ---------------------------------------------------------------------------- 
 
***************************************************************
* Done Running Example  ./ex2
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
