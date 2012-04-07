<?php $root=""; ?>
<?php require($root."navigation.php"); ?>
<html>
<head>
  <?php load_style($root); ?>
</head>
 
<body>
 
<?php make_navigation("reduced_basis_ex1",$root)?>
 
<div class="content">
<a name="comments"></a> 
<div class = "comment">
  

<br><br>rbOOmit is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.
  

<br><br>You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


<br><br>C++ include files that we need
</div>

<div class ="fragment">
<pre>
        #include &lt;iostream&gt;
        #include &lt;algorithm&gt;
        #include &lt;cmath&gt;
        #include &lt;set&gt;
        
</pre>
</div>
<div class = "comment">
Basic include file needed for the mesh functionality.
</div>

<div class ="fragment">
<pre>
        #include "libmesh.h"
        #include "mesh.h"
        #include "mesh_generation.h"
        #include "exodusII_io.h"
        #include "equation_systems.h"
        #include "dof_map.h"
        #include "getpot.h"
        #include "o_string_stream.h"
        #include "elem.h"
        
</pre>
</div>
<div class = "comment">
local includes
</div>

<div class ="fragment">
<pre>
        #include "rb_classes.h"
        #include "assembly.h"
        
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
<h1>Reduced Basis Example 1 - Certified Reduced Basis Method</h1>
 

<br><br>In this example problem we use the Certified Reduced Basis method
to solve a steady convection-diffusion problem on the unit square.
The reduced basis method relies on an expansion of the PDE in the
form \sum_q=1^Q_a theta_a^q(\mu) * a^q(u,v) = \sum_q=1^Q_f theta_f^q(\mu) f^q(v)
where theta_a, theta_f are parameter dependent functions and 
a^q, f^q are parameter independent operators (\mu denotes a parameter).


<br><br>We first attach the parameter dependent functions and paramater
independent operators to the RBSystem. Then in Offline mode, a
reduced basis space is generated and written out to the directory
"offline_data". In Online mode, the reduced basis data in "offline_data"
is read in and used to solve the reduced problem for the parameters
specified in reduced_basis_ex1.in.


<br><br>We also attach four outputs to the system which are averages over certain
subregions of the domain. In Online mode, we print out the values of these
outputs as well as rigorous error bounds with respect to the output
associated with the "truth" finite element discretization.


<br><br>The main program.
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
        
        #if !defined(LIBMESH_HAVE_XDR)
</pre>
</div>
<div class = "comment">
We need XDR support to write out reduced bases
</div>

<div class ="fragment">
<pre>
          libmesh_example_assert(false, "--enable-xdr");
        #elif defined(LIBMESH_DEFAULT_SINGLE_PRECISION)
</pre>
</div>
<div class = "comment">
XDR binary support requires double precision
</div>

<div class ="fragment">
<pre>
          libmesh_example_assert(false, "--disable-singleprecision");
        #endif
</pre>
</div>
<div class = "comment">
FIXME: This example currently segfaults with Trilinos?
</div>

<div class ="fragment">
<pre>
          libmesh_example_assert(libMesh::default_solver_package() == PETSC_SOLVERS, "--enable-petsc");
        
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
Parse the input file (reduced_basis_ex1.in) using GetPot
</div>

<div class ="fragment">
<pre>
          std::string parameters_filename = "reduced_basis_ex1.in";
          GetPot infile(parameters_filename);
        
          unsigned int n_elem = infile("n_elem", 1);       // Determines the number of elements in the "truth" mesh
          const unsigned int dim = 2;                      // The number of spatial dimensions
          
          bool store_basis_functions = infile("store_basis_functions", true); // Do we write the RB basis functions to disk?
        
</pre>
</div>
<div class = "comment">
Read the "online_mode" flag from the command line
</div>

<div class ="fragment">
<pre>
          GetPot command_line (argc, argv);
          int online_mode = 0;
          if ( command_line.search(1, "-online_mode") )
            online_mode = command_line.next(online_mode);
        
</pre>
</div>
<div class = "comment">
Build a mesh.
</div>

<div class ="fragment">
<pre>
          Mesh mesh (dim);
          MeshTools::Generation::build_square (mesh,
                                               n_elem, n_elem,
                                               0., 1.,
                                               0., 1.,
                                               QUAD4);
        
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
Build a new RBEvaluation object which will be used to perform
Reduced Basis calculations. This is required in both the
"Offline" and "Online" stages.
</div>

<div class ="fragment">
<pre>
          SimpleRBEvaluation rb_eval;
        
          if(!online_mode) // Perform the Offline stage of the RB method
          {
</pre>
</div>
<div class = "comment">
We override RBConstruction with SimpleRBConstruction in order to
specialize a few functions for this particular problem.
</div>

<div class ="fragment">
<pre>
            SimpleRBConstruction & rb_con =
              equation_systems.add_system&lt;SimpleRBConstruction&gt; ("RBConvectionDiffusion");
        
</pre>
</div>
<div class = "comment">
Initialize the data structures for the equation system.
</div>

<div class ="fragment">
<pre>
            equation_systems.init ();
        
</pre>
</div>
<div class = "comment">
Print out some information about the "truth" discretization
</div>

<div class ="fragment">
<pre>
            equation_systems.print_info();
            mesh.print_info();
          
</pre>
</div>
<div class = "comment">
We need to give the RBConstruction object a pointer to
our RBEvaluation object
</div>

<div class ="fragment">
<pre>
            rb_con.rb_eval = &rb_eval;
        
</pre>
</div>
<div class = "comment">
Read in the data that defines this problem from the specified text file
</div>

<div class ="fragment">
<pre>
            rb_con.process_parameters_file(parameters_filename);
        
</pre>
</div>
<div class = "comment">
Print out info that describes the current setup of rb_con
</div>

<div class ="fragment">
<pre>
            rb_con.print_info();
        
</pre>
</div>
<div class = "comment">
Prepare rb_con for the Construction stage of the RB method.
This sets up the necessary data structures and performs
initial assembly of the "truth" affine expansion of the PDE.
</div>

<div class ="fragment">
<pre>
            rb_con.initialize_rb_construction();
        
</pre>
</div>
<div class = "comment">
Compute the reduced basis space by computing "snapshots", i.e.
"truth" solves, at well-chosen parameter values and employing
these snapshots as basis functions.
</div>

<div class ="fragment">
<pre>
            rb_con.train_reduced_basis();
            
</pre>
</div>
<div class = "comment">
Write out the data that will subsequently be required for the Evaluation stage
</div>

<div class ="fragment">
<pre>
            rb_con.rb_eval-&gt;write_offline_data_to_files();
            
</pre>
</div>
<div class = "comment">
If requested, write out the RB basis functions for visualization purposes
</div>

<div class ="fragment">
<pre>
            if(store_basis_functions)
            {
</pre>
</div>
<div class = "comment">
If we want to be able to visualize the solution in the online stage,
then we should also save the state of the equation_systems object
so we can initialize it properly in the online stage
</div>

<div class ="fragment">
<pre>
              equation_systems.write("equation_systems.dat", WRITE);
              
</pre>
</div>
<div class = "comment">
Write out the basis functions
</div>

<div class ="fragment">
<pre>
              rb_con.rb_eval-&gt;write_out_basis_functions(rb_con);
            }
          }
          else // Perform the Online stage of the RB method
          {
</pre>
</div>
<div class = "comment">
Read in the reduced basis data
</div>

<div class ="fragment">
<pre>
            rb_eval.read_offline_data_from_files();
            
</pre>
</div>
<div class = "comment">
Get the parameters at which we do a reduced basis solve
</div>

<div class ="fragment">
<pre>
            unsigned int online_N = infile("online_N",1);
            unsigned int n_parameters = infile("n_parameters",1);
            std::vector&lt;Real&gt; online_mu_vector(n_parameters);
            for(unsigned int i=0; i&lt;n_parameters; i++)
            {
              online_mu_vector[i] = infile("online_mu", online_mu_vector[i], i);
            }
        
</pre>
</div>
<div class = "comment">
Set the parameters to online_mu_vector
</div>

<div class ="fragment">
<pre>
            rb_eval.set_current_parameters(online_mu_vector);
            rb_eval.print_current_parameters();
        
</pre>
</div>
<div class = "comment">
Now do the Online solve using the precomputed reduced basis
</div>

<div class ="fragment">
<pre>
            rb_eval.rb_solve(online_N);
        
</pre>
</div>
<div class = "comment">
Print out outputs as well as the corresponding output error bounds.
</div>

<div class ="fragment">
<pre>
            std::cout &lt;&lt; "output 1, value = " &lt;&lt; rb_eval.RB_outputs[0]
                      &lt;&lt; ", bound = " &lt;&lt; rb_eval.RB_output_error_bounds[0]
                      &lt;&lt; std::endl;
            std::cout &lt;&lt; "output 2, value = " &lt;&lt; rb_eval.RB_outputs[1]
                      &lt;&lt; ", bound = " &lt;&lt; rb_eval.RB_output_error_bounds[1]
                      &lt;&lt; std::endl;
            std::cout &lt;&lt; "output 3, value = " &lt;&lt; rb_eval.RB_outputs[2]
                      &lt;&lt; ", bound = " &lt;&lt; rb_eval.RB_output_error_bounds[2]
                      &lt;&lt; std::endl;
            std::cout &lt;&lt; "output 4, value = " &lt;&lt; rb_eval.RB_outputs[3]
                      &lt;&lt; ", bound = " &lt;&lt; rb_eval.RB_output_error_bounds[3]
                      &lt;&lt; std::endl &lt;&lt; std::endl;
        
            if(store_basis_functions)
            {
</pre>
</div>
<div class = "comment">
initialize the EquationSystems object by reading in the state that
was written out in the offline stage
</div>

<div class ="fragment">
<pre>
              equation_systems.read("equation_systems.dat", READ);
              RBConstruction& rb_con = equation_systems.get_system&lt;RBConstruction&gt;("RBConvectionDiffusion");
              rb_con.rb_eval = &rb_eval;
        
</pre>
</div>
<div class = "comment">
Read in the basis functions
</div>

<div class ="fragment">
<pre>
              rb_eval.read_in_basis_functions(rb_con);
              
</pre>
</div>
<div class = "comment">
Plot the solution
</div>

<div class ="fragment">
<pre>
              rb_con.load_rb_solution();
        #ifdef LIBMESH_HAVE_EXODUS_API
              ExodusII_IO(mesh).write_equation_systems ("RB_sol.e",equation_systems);
        #endif
              
</pre>
</div>
<div class = "comment">
Plot the first basis function that was generated from the train_reduced_basis
call in the Offline stage
</div>

<div class ="fragment">
<pre>
              rb_con.load_basis_function(0);
        #ifdef LIBMESH_HAVE_EXODUS_API
              ExodusII_IO(mesh).write_equation_systems ("bf0.e",equation_systems);
        #endif
            }
          }
        
          return 0;
        }
        
</pre>
</div>

<a name="nocomments"></a> 
<br><br><br> <h1> The program without comments: </h1> 
<pre> 
    
    
  
  #include &lt;iostream&gt;
  #include &lt;algorithm&gt;
  #include &lt;cmath&gt;
  #include &lt;set&gt;
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;mesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;mesh_generation.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;exodusII_io.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;equation_systems.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;dof_map.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;getpot.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;o_string_stream.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;elem.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;rb_classes.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;assembly.h&quot;</FONT></B>
  
  using namespace libMesh;
  
   
  
  
  
  <B><FONT COLOR="#228B22">int</FONT></B> main (<B><FONT COLOR="#228B22">int</FONT></B> argc, <B><FONT COLOR="#228B22">char</FONT></B>** argv)
  {
    LibMeshInit init (argc, argv);
  
  #<B><FONT COLOR="#A020F0">if</FONT></B> !defined(LIBMESH_HAVE_XDR)
    libmesh_example_assert(false, <B><FONT COLOR="#BC8F8F">&quot;--enable-xdr&quot;</FONT></B>);
  #elif defined(LIBMESH_DEFAULT_SINGLE_PRECISION)
    libmesh_example_assert(false, <B><FONT COLOR="#BC8F8F">&quot;--disable-singleprecision&quot;</FONT></B>);
  #endif
    libmesh_example_assert(libMesh::default_solver_package() == PETSC_SOLVERS, <B><FONT COLOR="#BC8F8F">&quot;--enable-petsc&quot;</FONT></B>);
  
    libmesh_example_assert(2 &lt;= LIBMESH_DIM, <B><FONT COLOR="#BC8F8F">&quot;2D support&quot;</FONT></B>);
    
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::string parameters_filename = <B><FONT COLOR="#BC8F8F">&quot;reduced_basis_ex1.in&quot;</FONT></B>;
    GetPot infile(parameters_filename);
  
    <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_elem = infile(<B><FONT COLOR="#BC8F8F">&quot;n_elem&quot;</FONT></B>, 1);       <I><FONT COLOR="#B22222">// Determines the number of elements in the &quot;truth&quot; mesh
</FONT></I>    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> dim = 2;                      <I><FONT COLOR="#B22222">// The number of spatial dimensions
</FONT></I>    
    <B><FONT COLOR="#228B22">bool</FONT></B> store_basis_functions = infile(<B><FONT COLOR="#BC8F8F">&quot;store_basis_functions&quot;</FONT></B>, true); <I><FONT COLOR="#B22222">// Do we write the RB basis functions to disk?
</FONT></I>  
    GetPot command_line (argc, argv);
    <B><FONT COLOR="#228B22">int</FONT></B> online_mode = 0;
    <B><FONT COLOR="#A020F0">if</FONT></B> ( command_line.search(1, <B><FONT COLOR="#BC8F8F">&quot;-online_mode&quot;</FONT></B>) )
      online_mode = command_line.next(online_mode);
  
    Mesh mesh (dim);
    <B><FONT COLOR="#5F9EA0">MeshTools</FONT></B>::Generation::build_square (mesh,
                                         n_elem, n_elem,
                                         0., 1.,
                                         0., 1.,
                                         QUAD4);
  
    EquationSystems equation_systems (mesh);
  
    SimpleRBEvaluation rb_eval;
  
    <B><FONT COLOR="#A020F0">if</FONT></B>(!online_mode) <I><FONT COLOR="#B22222">// Perform the Offline stage of the RB method
</FONT></I>    {
      SimpleRBConstruction &amp; rb_con =
        equation_systems.add_system&lt;SimpleRBConstruction&gt; (<B><FONT COLOR="#BC8F8F">&quot;RBConvectionDiffusion&quot;</FONT></B>);
  
      equation_systems.init ();
  
      equation_systems.print_info();
      mesh.print_info();
    
      rb_con.rb_eval = &amp;rb_eval;
  
      rb_con.process_parameters_file(parameters_filename);
  
      rb_con.print_info();
  
      rb_con.initialize_rb_construction();
  
      rb_con.train_reduced_basis();
      
      rb_con.rb_eval-&gt;write_offline_data_to_files();
      
      <B><FONT COLOR="#A020F0">if</FONT></B>(store_basis_functions)
      {
        equation_systems.write(<B><FONT COLOR="#BC8F8F">&quot;equation_systems.dat&quot;</FONT></B>, WRITE);
        
        rb_con.rb_eval-&gt;write_out_basis_functions(rb_con);
      }
    }
    <B><FONT COLOR="#A020F0">else</FONT></B> <I><FONT COLOR="#B22222">// Perform the Online stage of the RB method
</FONT></I>    {
      rb_eval.read_offline_data_from_files();
      
      <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> online_N = infile(<B><FONT COLOR="#BC8F8F">&quot;online_N&quot;</FONT></B>,1);
      <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_parameters = infile(<B><FONT COLOR="#BC8F8F">&quot;n_parameters&quot;</FONT></B>,1);
      <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;Real&gt; online_mu_vector(n_parameters);
      <B><FONT COLOR="#A020F0">for</FONT></B>(<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i&lt;n_parameters; i++)
      {
        online_mu_vector[i] = infile(<B><FONT COLOR="#BC8F8F">&quot;online_mu&quot;</FONT></B>, online_mu_vector[i], i);
      }
  
      rb_eval.set_current_parameters(online_mu_vector);
      rb_eval.print_current_parameters();
  
      rb_eval.rb_solve(online_N);
  
      <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;output 1, value = &quot;</FONT></B> &lt;&lt; rb_eval.RB_outputs[0]
                &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;, bound = &quot;</FONT></B> &lt;&lt; rb_eval.RB_output_error_bounds[0]
                &lt;&lt; std::endl;
      <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;output 2, value = &quot;</FONT></B> &lt;&lt; rb_eval.RB_outputs[1]
                &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;, bound = &quot;</FONT></B> &lt;&lt; rb_eval.RB_output_error_bounds[1]
                &lt;&lt; std::endl;
      <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;output 3, value = &quot;</FONT></B> &lt;&lt; rb_eval.RB_outputs[2]
                &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;, bound = &quot;</FONT></B> &lt;&lt; rb_eval.RB_output_error_bounds[2]
                &lt;&lt; std::endl;
      <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;output 4, value = &quot;</FONT></B> &lt;&lt; rb_eval.RB_outputs[3]
                &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;, bound = &quot;</FONT></B> &lt;&lt; rb_eval.RB_output_error_bounds[3]
                &lt;&lt; std::endl &lt;&lt; std::endl;
  
      <B><FONT COLOR="#A020F0">if</FONT></B>(store_basis_functions)
      {
        equation_systems.read(<B><FONT COLOR="#BC8F8F">&quot;equation_systems.dat&quot;</FONT></B>, READ);
        RBConstruction&amp; rb_con = equation_systems.get_system&lt;RBConstruction&gt;(<B><FONT COLOR="#BC8F8F">&quot;RBConvectionDiffusion&quot;</FONT></B>);
        rb_con.rb_eval = &amp;rb_eval;
  
        rb_eval.read_in_basis_functions(rb_con);
        
        rb_con.load_rb_solution();
  #ifdef LIBMESH_HAVE_EXODUS_API
        ExodusII_IO(mesh).write_equation_systems (<B><FONT COLOR="#BC8F8F">&quot;RB_sol.e&quot;</FONT></B>,equation_systems);
  #endif
        
        rb_con.load_basis_function(0);
  #ifdef LIBMESH_HAVE_EXODUS_API
        ExodusII_IO(mesh).write_equation_systems (<B><FONT COLOR="#BC8F8F">&quot;bf0.e&quot;</FONT></B>,equation_systems);
  #endif
      }
    }
  
    <B><FONT COLOR="#A020F0">return</FONT></B> 0;
  }
  
</pre> 
<a name="output"></a> 
<br><br><br> <h1> The console output of the program: </h1> 
<pre>
Updated .depend
Compiling C++ (in optimized mode) reduced_basis_ex1.C...
Linking reduced_basis_ex1-opt...
***************************************************************
* Running  ./reduced_basis_ex1-opt
***************************************************************
 
*** Warning, This code is untested, experimental, or likely to see future API changes: src/reduced_basis/rb_parametrized_object.C, line 31, compiled Apr  7 2012 at 15:51:42 ***
 EquationSystems
  n_systems()=1
   System #0, "RBConvectionDiffusion"
    Type "RBConstruction"
    Variables="u" 
    Finite Element Types="LAGRANGE", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" 
    Approximation Orders="FIRST", "THIRD" 
    n_dofs()=2601
    n_local_dofs()=2601
    n_constrained_dofs()=200
    n_local_constrained_dofs()=200
    n_vectors()=1
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 8.76624
      Average Off-Processor Bandwidth <= 0
      Maximum  On-Processor Bandwidth <= 9
      Maximum Off-Processor Bandwidth <= 0
    DofMap Constraints
      Number of DoF Constraints = 200
      Average DoF Constraint Length= 0
      Number of Node Constraints = 0

 Mesh Information:
  mesh_dimension()=2
  spatial_dimension()=3
  n_nodes()=2601
    n_local_nodes()=2601
  n_elem()=2500
    n_local_elem()=2500
    n_active_elem()=2500
  n_subdomains()=1
  n_partitions()=1
  n_processors()=1
  n_threads()=1
  processor_id()=0

Initializing training parameters with deterministic training set...
Parameter 0: log scaling = 0
Parameter 1: log scaling = 0


RBConstruction parameters:
system name: RBConvectionDiffusion
constrained_problem: 0
Nmax: 20
A_q operators attached: 3
F_q functions attached: 1
n_outputs: 4
output 0, Q_l = 1
output 1, Q_l = 1
output 2, Q_l = 1
output 3, Q_l = 1
Number of parameters: 2
Parameter 0: Min = -2, Max = 2
Parameter 1: Min = -2, Max = 2
n_training_samples: 100
single-matrix mode? 0
reuse preconditioner? 1
use a relative error bound in greedy? 0
write out data during basis training? 0
quiet mode? 1
parameter initialized to: 
mu[0] = -2
mu[1] = -2

Compute output dual norms
output_dual_norms[0][0] = 0.323675
output_dual_norms[1][0] = 0.323675
output_dual_norms[2][0] = 0.323675
output_dual_norms[3][0] = 0.323675

---- Performing Greedy basis enrichment ----

---- Basis dimension: 0 ----
Performing RB solves on training set
Maximum (absolute) error bound is 3.74824

Performing truth solve at parameter:
mu[0] = -2
mu[1] = -2

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 1 ----
Performing RB solves on training set
Maximum (absolute) error bound is 6.67625

Performing truth solve at parameter:
mu[0] = 2
mu[1] = 2

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 2 ----
Performing RB solves on training set
Maximum (absolute) error bound is 5.18295

Performing truth solve at parameter:
mu[0] = -2
mu[1] = 2

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 3 ----
Performing RB solves on training set
Maximum (absolute) error bound is 3.59028

Performing truth solve at parameter:
mu[0] = 2
mu[1] = -2

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 4 ----
Performing RB solves on training set
Maximum (absolute) error bound is 2.71338

Performing truth solve at parameter:
mu[0] = 0.222222
mu[1] = 0.222222

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 5 ----
Performing RB solves on training set
Maximum (absolute) error bound is 2.50784

Performing truth solve at parameter:
mu[0] = -0.666667
mu[1] = -0.222222

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 6 ----
Performing RB solves on training set
Maximum (absolute) error bound is 2.15653

Performing truth solve at parameter:
mu[0] = 0.222222
mu[1] = -0.666667

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 7 ----
Performing RB solves on training set
Maximum (absolute) error bound is 1.71352

Performing truth solve at parameter:
mu[0] = -0.222222
mu[1] = 1.11111

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 8 ----
Performing RB solves on training set
Maximum (absolute) error bound is 1.61442

Performing truth solve at parameter:
mu[0] = 1.55556
mu[1] = -0.222222

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 9 ----
Performing RB solves on training set
Maximum (absolute) error bound is 1.09712

Performing truth solve at parameter:
mu[0] = -2
mu[1] = 0.222222

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 10 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.9794

Performing truth solve at parameter:
mu[0] = -0.222222
mu[1] = -2

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 11 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.763958

Performing truth solve at parameter:
mu[0] = -0.222222
mu[1] = 0.222222

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 12 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.616746

Performing truth solve at parameter:
mu[0] = 1.11111
mu[1] = 1.11111

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 13 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.530181

Performing truth solve at parameter:
mu[0] = 0.222222
mu[1] = 2

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 14 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.475932

Performing truth solve at parameter:
mu[0] = 0.666667
mu[1] = -0.222222

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 15 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.440399

Performing truth solve at parameter:
mu[0] = -1.11111
mu[1] = 0.666667

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 16 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.423537

Performing truth solve at parameter:
mu[0] = 1.11111
mu[1] = -1.11111

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 17 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.408558

Performing truth solve at parameter:
mu[0] = -0.666667
mu[1] = -1.11111

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 18 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.333764

Performing truth solve at parameter:
mu[0] = -0.222222
mu[1] = -0.222222

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 19 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.196867

Performing truth solve at parameter:
mu[0] = 0.222222
mu[1] = 0.666667

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 20 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.191475

Maximum number of basis functions reached: Nmax = 20
Writing out the basis functions...

-------------------------------------------------------------------
| Time:           Sat Apr  7 16:00:57 2012                         |
| OS:             Linux                                            |
| HostName:       lkirk-home                                       |
| OS Release:     3.0.0-17-generic                                 |
| OS Version:     #30-Ubuntu SMP Thu Mar 8 20:45:39 UTC 2012       |
| Machine:        x86_64                                           |
| Username:       benkirk                                          |
| Configuration:  ./configure run on Sat Apr  7 15:49:27 CDT 2012  |
-------------------------------------------------------------------
 --------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=3.32218, Active time=3.01249                                                 |
 --------------------------------------------------------------------------------------------------------------
| Event                            nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                            w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|--------------------------------------------------------------------------------------------------------------|
|                                                                                                              |
|                                                                                                              |
| DofMap                                                                                                       |
|   add_neighbors_to_send_list()   1         0.0014      0.001403    0.0014      0.001403    0.05     0.05     |
|   build_constraint_matrix()      22500     0.0138      0.000001    0.0138      0.000001    0.46     0.46     |
|   cnstrn_elem_mat_vec()          22500     0.0295      0.000001    0.0295      0.000001    0.98     0.98     |
|   compute_sparsity()             1         0.0071      0.007074    0.0131      0.013052    0.23     0.43     |
|   create_dof_constraints()       1         0.0086      0.008613    0.0121      0.012079    0.29     0.40     |
|   distribute_dofs()              1         0.0021      0.002113    0.0083      0.008334    0.07     0.28     |
|   dof_indices()                  50000     0.0476      0.000001    0.0476      0.000001    1.58     1.58     |
|   prepare_send_list()            1         0.0000      0.000001    0.0000      0.000001    0.00     0.00     |
|   reinit()                       1         0.0062      0.006218    0.0062      0.006218    0.21     0.21     |
|                                                                                                              |
| EquationSystems                                                                                              |
|   write()                        1         0.0147      0.014702    0.0151      0.015058    0.49     0.50     |
|                                                                                                              |
| FE                                                                                                           |
|   compute_affine_map()           24300     0.0180      0.000001    0.0180      0.000001    0.60     0.60     |
|   compute_face_map()             1800      0.0055      0.000003    0.0118      0.000007    0.18     0.39     |
|   compute_shape_functions()      24300     0.0117      0.000000    0.0117      0.000000    0.39     0.39     |
|   init_face_shape_functions()    9         0.0000      0.000003    0.0000      0.000003    0.00     0.00     |
|   init_shape_functions()         1809      0.0078      0.000004    0.0078      0.000004    0.26     0.26     |
|   inverse_map()                  3600      0.0056      0.000002    0.0056      0.000002    0.19     0.19     |
|                                                                                                              |
| Mesh                                                                                                         |
|   find_neighbors()               1         0.0078      0.007810    0.0078      0.007810    0.26     0.26     |
|   renumber_nodes_and_elem()      2         0.0007      0.000326    0.0007      0.000326    0.02     0.02     |
|                                                                                                              |
| MeshCommunication                                                                                            |
|   assign_global_indices()        1         0.0509      0.050876    0.0509      0.050910    1.69     1.69     |
|                                                                                                              |
| MeshTools::Generation                                                                                        |
|   build_cube()                   1         0.0041      0.004102    0.0041      0.004102    0.14     0.14     |
|                                                                                                              |
| Parallel                                                                                                     |
|   allgather()                    5         0.0000      0.000000    0.0000      0.000000    0.00     0.00     |
|   receive()                      84        0.0003      0.000004    0.0003      0.000004    0.01     0.01     |
|   send()                         84        0.0009      0.000011    0.0009      0.000011    0.03     0.03     |
|   send_receive()                 4         0.0000      0.000007    0.0000      0.000007    0.00     0.00     |
|                                                                                                              |
| Partitioner                                                                                                  |
|   single_partition()             1         0.0004      0.000372    0.0004      0.000372    0.01     0.01     |
|                                                                                                              |
| PetscLinearSolver                                                                                            |
|   solve()                        85        1.7594      0.020699    1.7594      0.020699    58.40    58.40    |
|                                                                                                              |
| RBConstruction                                                                                               |
|   add_scaled_matrix_and_vector() 9         0.1615      0.017948    0.3166      0.035173    5.36     10.51    |
|   clear()                        1         0.0004      0.000403    0.0004      0.000403    0.01     0.01     |
|   compute_Fq_representor_norms() 1         0.0011      0.001054    0.0202      0.020173    0.03     0.67     |
|   compute_max_error_bound()      21        0.0097      0.000462    0.0971      0.004624    0.32     3.22     |
|   compute_output_dual_norms()    1         0.1690      0.168958    0.2610      0.260952    5.61     8.66     |
|   enrich_RB_space()              20        0.0298      0.001489    0.0298      0.001489    0.99     0.99     |
|   train_reduced_basis()          1         0.0038      0.003792    2.6010      2.600975    0.13     86.34    |
|   truth_assembly()               20        0.0521      0.002603    0.0521      0.002603    1.73     1.73     |
|   truth_solve()                  20        0.0052      0.000259    0.3732      0.018658    0.17     12.39    |
|   update_RB_system_matrices()    20        0.1570      0.007850    0.1570      0.007850    5.21     5.21     |
|   update_residual_terms()        20        0.3265      0.016324    1.6589      0.082943    10.84    55.07    |
|                                                                                                              |
| RBEvaluation                                                                                                 |
|   clear()                        1         0.0001      0.000080    0.0002      0.000210    0.00     0.01     |
|   clear_riesz_representors()     2         0.0001      0.000064    0.0001      0.000064    0.00     0.00     |
|   compute_residual_dual_norm()   2100      0.0576      0.000027    0.0576      0.000027    1.91     1.91     |
|   rb_solve()                     2100      0.0291      0.000014    0.0871      0.000041    0.97     2.89     |
|   resize_data_structures()       1         0.0001      0.000083    0.0001      0.000083    0.00     0.00     |
|   write_offline_data_to_files()  1         0.0055      0.005486    0.0055      0.005486    0.18     0.18     |
 --------------------------------------------------------------------------------------------------------------
| Totals:                          155432    3.0125                                          100.00            |
 --------------------------------------------------------------------------------------------------------------

*** Warning, This code is untested, experimental, or likely to see future API changes: src/reduced_basis/rb_parametrized_object.C, line 31, compiled Apr  7 2012 at 15:51:42 ***
mu[0] = -2
mu[1] = 0.79

output 1, value = 0.12248, bound = 0.0659666
output 2, value = 0.370889, bound = 0.0659666
output 3, value = 0.254934, bound = 0.0659666
output 4, value = 0.114764, bound = 0.0659666

Reading in the basis functions...
Finished reading in the basis functions...

-------------------------------------------------------------------
| Time:           Sat Apr  7 16:00:58 2012                         |
| OS:             Linux                                            |
| HostName:       lkirk-home                                       |
| OS Release:     3.0.0-17-generic                                 |
| OS Version:     #30-Ubuntu SMP Thu Mar 8 20:45:39 UTC 2012       |
| Machine:        x86_64                                           |
| Username:       benkirk                                          |
| Configuration:  ./configure run on Sat Apr  7 15:49:27 CDT 2012  |
-------------------------------------------------------------------
 --------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=0.367923, Active time=0.128905                                               |
 --------------------------------------------------------------------------------------------------------------
| Event                            nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                            w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|--------------------------------------------------------------------------------------------------------------|
|                                                                                                              |
|                                                                                                              |
| DofMap                                                                                                       |
|   add_neighbors_to_send_list()   1         0.0007      0.000702    0.0007      0.000702    0.54     0.54     |
|   compute_sparsity()             1         0.0052      0.005184    0.0095      0.009525    4.02     7.39     |
|   create_dof_constraints()       1         0.0006      0.000611    0.0006      0.000611    0.47     0.47     |
|   distribute_dofs()              1         0.0012      0.001154    0.0060      0.005954    0.90     4.62     |
|   dof_indices()                  7500      0.0087      0.000001    0.0087      0.000001    6.71     6.71     |
|   prepare_send_list()            1         0.0000      0.000001    0.0000      0.000001    0.00     0.00     |
|   reinit()                       1         0.0048      0.004799    0.0048      0.004799    3.72     3.72     |
|                                                                                                              |
| EquationSystems                                                                                              |
|   build_solution_vector()        2         0.0064      0.003206    0.0120      0.006000    4.97     9.31     |
|   read()                         1         0.0160      0.015952    0.0892      0.089240    12.38    69.23    |
|   update()                       1         0.0001      0.000062    0.0001      0.000062    0.05     0.05     |
|                                                                                                              |
| ExodusII_IO                                                                                                  |
|   write_nodal_data()             2         0.0069      0.003429    0.0069      0.003429    5.32     5.32     |
|                                                                                                              |
| Mesh                                                                                                         |
|   find_neighbors()               1         0.0077      0.007681    0.0077      0.007681    5.96     5.96     |
|   renumber_nodes_and_elem()      2         0.0006      0.000325    0.0006      0.000325    0.50     0.50     |
|                                                                                                              |
| MeshCommunication                                                                                            |
|   assign_global_indices()        1         0.0564      0.056407    0.0564      0.056433    43.76    43.78    |
|                                                                                                              |
| MeshOutput                                                                                                   |
|   write_equation_systems()       2         0.0000      0.000025    0.0189      0.009454    0.04     14.67    |
|                                                                                                              |
| MeshTools::Generation                                                                                        |
|   build_cube()                   1         0.0040      0.004009    0.0040      0.004009    3.11     3.11     |
|                                                                                                              |
| Parallel                                                                                                     |
|   allgather()                    47        0.0001      0.000001    0.0001      0.000001    0.05     0.05     |
|   send_receive()                 4         0.0000      0.000005    0.0000      0.000005    0.02     0.02     |
|                                                                                                              |
| Partitioner                                                                                                  |
|   single_partition()             1         0.0004      0.000352    0.0004      0.000352    0.27     0.27     |
|                                                                                                              |
| RBConstruction                                                                                               |
|   clear()                        2         0.0002      0.000092    0.0002      0.000092    0.14     0.14     |
|   load_basis_function()          1         0.0000      0.000045    0.0000      0.000045    0.03     0.03     |
|   load_rb_solution()             1         0.0002      0.000205    0.0002      0.000205    0.16     0.16     |
|                                                                                                              |
| RBEvaluation                                                                                                 |
|   clear()                        1         0.0001      0.000144    0.0001      0.000146    0.11     0.11     |
|   clear_riesz_representors()     2         0.0000      0.000003    0.0000      0.000003    0.00     0.00     |
|   compute_residual_dual_norm()   1         0.0002      0.000177    0.0002      0.000177    0.14     0.14     |
|   rb_solve()                     1         0.0006      0.000581    0.0008      0.000758    0.45     0.59     |
|   read_offline_data_from_files() 1         0.0078      0.007759    0.0080      0.007951    6.02     6.17     |
|   resize_data_structures()       1         0.0002      0.000186    0.0002      0.000192    0.14     0.15     |
 --------------------------------------------------------------------------------------------------------------
| Totals:                          7582      0.1289                                          100.00            |
 --------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running  ./reduced_basis_ex1-opt
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
