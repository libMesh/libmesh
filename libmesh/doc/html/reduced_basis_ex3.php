<?php $root=""; ?>
<?php require($root."navigation.php"); ?>
<html>
<head>
  <?php load_style($root); ?>
</head>
 
<body>
 
<?php make_navigation("reduced_basis_ex3",$root)?>
 
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


<br><br>Basic include file needed for the mesh functionality.
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
<h1>Reduced Basis Example 3 - Transient Reduced Basis Problem</h1>


<br><br>In this example problem we use the Certified Reduced Basis method
to solve a transient convection-diffusion problem on the unit square.
The PDE is similar to reduced_basis_ex1, except there is a time-derivative
in this case.


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
        
</pre>
</div>
<div class = "comment">
This example requires SLEPc
</div>

<div class ="fragment">
<pre>
        #if !defined(LIBMESH_HAVE_SLEPC)
          libmesh_example_assert(false, "--enable-slepc");
        #else
        
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
Parse the input file (reduced_basis_ex3.in) using GetPot
</div>

<div class ="fragment">
<pre>
          std::string parameters_filename = "reduced_basis_ex3.in";
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
Finally, we need to give the RBConstruction object a pointer to
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
            rb_eval.temporal_discretization = rb_con.temporal_discretization;
        
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
            rb_eval.process_temporal_parameters_file(parameters_filename);
            rb_eval.print_current_parameters();
        
</pre>
</div>
<div class = "comment">
Now do the Online solve using the precomputed reduced basis
</div>

<div class ="fragment">
<pre>
            Real error_bound_final_time = rb_eval.rb_solve(online_N);
            
            libMesh::out &lt;&lt; "Error bound (absolute) at the final time is "
                         &lt;&lt; error_bound_final_time &lt;&lt; std::endl &lt;&lt; std::endl;
        
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
              TransientRBConstruction& rb_con =
                equation_systems.get_system&lt;TransientRBConstruction&gt;("RBConvectionDiffusion");
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
Plot the solution at the final time level
</div>

<div class ="fragment">
<pre>
              rb_con.temporal_discretization = rb_eval.temporal_discretization;
              rb_con.temporal_discretization.set_time_step(rb_con.temporal_discretization.get_n_time_steps());
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
        
        #endif // LIBMESH_HAVE_SLEPC
        }
        
</pre>
</div>

<a name="nocomments"></a> 
<br><br><br> <h1> The program without comments: </h1> 
<pre> 
    
    
  
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
  
  #<B><FONT COLOR="#A020F0">if</FONT></B> !defined(LIBMESH_HAVE_SLEPC)
    libmesh_example_assert(false, <B><FONT COLOR="#BC8F8F">&quot;--enable-slepc&quot;</FONT></B>);
  #<B><FONT COLOR="#A020F0">else</FONT></B>
  
  #<B><FONT COLOR="#A020F0">if</FONT></B> !defined(LIBMESH_HAVE_XDR)
    libmesh_example_assert(false, <B><FONT COLOR="#BC8F8F">&quot;--enable-xdr&quot;</FONT></B>);
  #elif defined(LIBMESH_DEFAULT_SINGLE_PRECISION)
    libmesh_example_assert(false, <B><FONT COLOR="#BC8F8F">&quot;--disable-singleprecision&quot;</FONT></B>);
  #endif
    libmesh_example_assert(libMesh::default_solver_package() == PETSC_SOLVERS, <B><FONT COLOR="#BC8F8F">&quot;--enable-petsc&quot;</FONT></B>);
  
    libmesh_example_assert(2 &lt;= LIBMESH_DIM, <B><FONT COLOR="#BC8F8F">&quot;2D support&quot;</FONT></B>);
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::string parameters_filename = <B><FONT COLOR="#BC8F8F">&quot;reduced_basis_ex3.in&quot;</FONT></B>;
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
      rb_eval.temporal_discretization = rb_con.temporal_discretization;
  
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
      rb_eval.process_temporal_parameters_file(parameters_filename);
      rb_eval.print_current_parameters();
  
      Real error_bound_final_time = rb_eval.rb_solve(online_N);
      
      <B><FONT COLOR="#5F9EA0">libMesh</FONT></B>::out &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Error bound (absolute) at the final time is &quot;</FONT></B>
                   &lt;&lt; error_bound_final_time &lt;&lt; std::endl &lt;&lt; std::endl;
  
      <B><FONT COLOR="#A020F0">if</FONT></B>(store_basis_functions)
      {
        equation_systems.read(<B><FONT COLOR="#BC8F8F">&quot;equation_systems.dat&quot;</FONT></B>, READ);
        TransientRBConstruction&amp; rb_con =
          equation_systems.get_system&lt;TransientRBConstruction&gt;(<B><FONT COLOR="#BC8F8F">&quot;RBConvectionDiffusion&quot;</FONT></B>);
        rb_con.rb_eval = &amp;rb_eval;
  
        rb_eval.read_in_basis_functions(rb_con);
        
        rb_con.temporal_discretization = rb_eval.temporal_discretization;
        rb_con.temporal_discretization.set_time_step(rb_con.temporal_discretization.get_n_time_steps());
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
  
  #endif <I><FONT COLOR="#B22222">// LIBMESH_HAVE_SLEPC
</FONT></I>  }
  
</pre> 
<a name="output"></a> 
<br><br><br> <h1> The console output of the program: </h1> 
<pre>
Updated .depend
Compiling C++ (in optimized mode) reduced_basis_ex3.C...
Linking reduced_basis_ex3-opt...
***************************************************************
* Running  ./reduced_basis_ex3-opt
***************************************************************
 
*** Warning, This code is untested, experimental, or likely to see future API changes: src/reduced_basis/rb_parametrized_object.C, line 31, compiled Apr  7 2012 at 15:51:42 ***
 EquationSystems
  n_systems()=1
   System #0, "RBConvectionDiffusion"
    Type "TransientRBConstruction"
    Variables="u" 
    Finite Element Types="LAGRANGE", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" 
    Approximation Orders="FIRST", "THIRD" 
    n_dofs()=676
    n_local_dofs()=676
    n_constrained_dofs()=100
    n_local_constrained_dofs()=100
    n_vectors()=3
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 8.54438
      Average Off-Processor Bandwidth <= 0
      Maximum  On-Processor Bandwidth <= 9
      Maximum Off-Processor Bandwidth <= 0
    DofMap Constraints
      Number of DoF Constraints = 100
      Average DoF Constraint Length= 0
      Number of Node Constraints = 0

 Mesh Information:
  mesh_dimension()=2
  spatial_dimension()=3
  n_nodes()=676
    n_local_nodes()=676
  n_elem()=625
    n_local_elem()=625
    n_active_elem()=625
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
Parameter 0: Min = 0, Max = 2
Parameter 1: Min = 0, Max = 2
n_training_samples: 100
single-matrix mode? 0
reuse preconditioner? 1
use a relative error bound in greedy? 1
write out data during basis training? 0
quiet mode? 1
parameter initialized to: 
mu[0] = 0
mu[1] = 0


TransientRBConstruction parameters:
Q_m: 1
Number of time-steps: 100
dt: 0.01
euler_theta (time discretization parameter): 1
delta_N (number of basis functions to add each POD-Greedy step): 1
Using zero initial condition

Compute output dual norms
output_dual_norms[0][0] = 33.6538
output_dual_norms[1][0] = 33.6538
output_dual_norms[2][0] = 33.6538
output_dual_norms[3][0] = 33.6538

---- Performing Greedy basis enrichment ----

---- Basis dimension: 0 ----
Performing RB solves on training set
Maximum (relative) error bound is inf

Performing truth solve at parameter:
mu[0] = 0
mu[1] = 0

Enriching the RB space

POD Eigenvalues:
eigenvalue 0 = 2.69179
eigenvalue 1 = 0.0276413
eigenvalue 2 = 0.000875017
...
last eigenvalue = -1.13824e-16

Updating RB matrices
Updating RB residual terms
Updating RB initial conditions

---- Basis dimension: 1 ----
Performing RB solves on training set
Maximum (relative) error bound is 2.58088

Performing truth solve at parameter:
mu[0] = 2
mu[1] = 2

Enriching the RB space

POD Eigenvalues:
eigenvalue 0 = 1.55804
eigenvalue 1 = 0.00755328
eigenvalue 2 = 0.00103931
...
last eigenvalue = -3.21128e-16

Updating RB matrices
Updating RB residual terms
Updating RB initial conditions

---- Basis dimension: 2 ----
Performing RB solves on training set
Maximum (relative) error bound is 1.46486

Performing truth solve at parameter:
mu[0] = 2
mu[1] = 0

Enriching the RB space

POD Eigenvalues:
eigenvalue 0 = 0.952676
eigenvalue 1 = 0.00511887
eigenvalue 2 = 0.00103315
...
last eigenvalue = -4.3702e-17

Updating RB matrices
Updating RB residual terms
Updating RB initial conditions

---- Basis dimension: 3 ----
Performing RB solves on training set
Maximum (relative) error bound is 1.47559

Performing truth solve at parameter:
mu[0] = 0
mu[1] = 2

Enriching the RB space

POD Eigenvalues:
eigenvalue 0 = 0.659245
eigenvalue 1 = 0.00491507
eigenvalue 2 = 0.000974904
...
last eigenvalue = -2.37661e-17

Updating RB matrices
Updating RB residual terms
Updating RB initial conditions

---- Basis dimension: 4 ----
Performing RB solves on training set
Maximum (relative) error bound is 0.516964

Performing truth solve at parameter:
mu[0] = 0.888889
mu[1] = 1.11111

Enriching the RB space

POD Eigenvalues:
eigenvalue 0 = 0.141459
eigenvalue 1 = 0.00444052
eigenvalue 2 = 0.00204462
...
last eigenvalue = -1.87763e-17

Updating RB matrices
Updating RB residual terms
Updating RB initial conditions

---- Basis dimension: 5 ----
Performing RB solves on training set
Maximum (relative) error bound is 0.302528

Performing truth solve at parameter:
mu[0] = 1.11111
mu[1] = 0

Enriching the RB space

POD Eigenvalues:
eigenvalue 0 = 0.0639292
eigenvalue 1 = 0.007257
eigenvalue 2 = 0.00169663
...
last eigenvalue = -1.59913e-18

Updating RB matrices
Updating RB residual terms
Updating RB initial conditions

---- Basis dimension: 6 ----
Performing RB solves on training set
Maximum (relative) error bound is 0.288478

Performing truth solve at parameter:
mu[0] = 0.888889
mu[1] = 2

Enriching the RB space

POD Eigenvalues:
eigenvalue 0 = 0.0196529
eigenvalue 1 = 0.00225713
eigenvalue 2 = 0.000843918
...
last eigenvalue = -4.51223e-19

Updating RB matrices
Updating RB residual terms
Updating RB initial conditions

---- Basis dimension: 7 ----
Performing RB solves on training set
Maximum (relative) error bound is 0.267677

Performing truth solve at parameter:
mu[0] = 2
mu[1] = 0.888889

Enriching the RB space

POD Eigenvalues:
eigenvalue 0 = 0.0108148
eigenvalue 1 = 0.00208092
eigenvalue 2 = 0.000703062
...
last eigenvalue = -1.22406e-18

Updating RB matrices
Updating RB residual terms
Updating RB initial conditions

---- Basis dimension: 8 ----
Performing RB solves on training set
Maximum (relative) error bound is 0.203171

Performing truth solve at parameter:
mu[0] = 0
mu[1] = 0.888889

Enriching the RB space

POD Eigenvalues:
eigenvalue 0 = 0.0410599
eigenvalue 1 = 0.006897
eigenvalue 2 = 0.00147376
...
last eigenvalue = -6.03973e-19

Updating RB matrices
Updating RB residual terms
Updating RB initial conditions

---- Basis dimension: 9 ----
Performing RB solves on training set
Maximum (relative) error bound is 0.121711

Performing truth solve at parameter:
mu[0] = 2
mu[1] = 1.55556

Enriching the RB space

POD Eigenvalues:
eigenvalue 0 = 0.00157141
eigenvalue 1 = 0.000569085
eigenvalue 2 = 0.000297784
...
last eigenvalue = -4.5915e-20

Updating RB matrices
Updating RB residual terms
Updating RB initial conditions

---- Basis dimension: 10 ----
Performing RB solves on training set
Maximum (relative) error bound is 0.121894

Performing truth solve at parameter:
mu[0] = 0
mu[1] = 1.55556

Enriching the RB space

POD Eigenvalues:
eigenvalue 0 = 0.003803
eigenvalue 1 = 0.00220747
eigenvalue 2 = 0.000457253
...
last eigenvalue = -3.25261e-19

Updating RB matrices
Updating RB residual terms
Updating RB initial conditions

---- Basis dimension: 11 ----
Performing RB solves on training set
Maximum (relative) error bound is 0.115558

Performing truth solve at parameter:
mu[0] = 1.33333
mu[1] = 2

Enriching the RB space

POD Eigenvalues:
eigenvalue 0 = 0.000968525
eigenvalue 1 = 0.000341279
eigenvalue 2 = 0.000319966
...
last eigenvalue = -2.11406e-20

Updating RB matrices
Updating RB residual terms
Updating RB initial conditions

---- Basis dimension: 12 ----
Performing RB solves on training set
Maximum (relative) error bound is 0.0901998

Performing truth solve at parameter:
mu[0] = 1.55556
mu[1] = 0.444444

Enriching the RB space

POD Eigenvalues:
eigenvalue 0 = 0.00247619
eigenvalue 1 = 0.00137033
eigenvalue 2 = 0.000358375
...
last eigenvalue = -1.69407e-19

Updating RB matrices
Updating RB residual terms
Updating RB initial conditions

---- Basis dimension: 13 ----
Performing RB solves on training set
Maximum (relative) error bound is 0.0616133

Performing truth solve at parameter:
mu[0] = 1.55556
mu[1] = 1.55556

Enriching the RB space

POD Eigenvalues:
eigenvalue 0 = 0.000524831
eigenvalue 1 = 0.000175057
eigenvalue 2 = 7.21479e-05
...
last eigenvalue = -8.45625e-21

Updating RB matrices
Updating RB residual terms
Updating RB initial conditions

---- Basis dimension: 14 ----
Performing RB solves on training set
Maximum (relative) error bound is 0.0578913

Performing truth solve at parameter:
mu[0] = 1.55556
mu[1] = 0.444444

Enriching the RB space

POD Eigenvalues:
eigenvalue 0 = 0.00136485
eigenvalue 1 = 0.000357777
eigenvalue 2 = 9.30974e-05
...
last eigenvalue = -4.74338e-20

Updating RB matrices
Updating RB residual terms
Updating RB initial conditions

---- Basis dimension: 15 ----
Performing RB solves on training set
Maximum (relative) error bound is 0.0511093

Performing truth solve at parameter:
mu[0] = 0
mu[1] = 1.55556

Enriching the RB space

POD Eigenvalues:
eigenvalue 0 = 0.000890358
eigenvalue 1 = 0.000289399
eigenvalue 2 = 8.7165e-05
...
last eigenvalue = -2.66207e-20

Updating RB matrices
Updating RB residual terms
Updating RB initial conditions

---- Basis dimension: 16 ----
Performing RB solves on training set
Maximum (relative) error bound is 0.0464637

Performing truth solve at parameter:
mu[0] = 2
mu[1] = 0

Enriching the RB space

POD Eigenvalues:
eigenvalue 0 = 0.000324742
eigenvalue 1 = 0.000211785
eigenvalue 2 = 4.51171e-05
...
last eigenvalue = -9.97465e-21

Updating RB matrices
Updating RB residual terms
Updating RB initial conditions

---- Basis dimension: 17 ----
Performing RB solves on training set
Maximum (relative) error bound is 0.0361278

Performing truth solve at parameter:
mu[0] = 0.222222
mu[1] = 2

Enriching the RB space

POD Eigenvalues:
eigenvalue 0 = 0.000312122
eigenvalue 1 = 7.29093e-05
eigenvalue 2 = 5.59736e-05
...
last eigenvalue = -8.73865e-21

Updating RB matrices
Updating RB residual terms
Updating RB initial conditions

---- Basis dimension: 18 ----
Performing RB solves on training set
Maximum (relative) error bound is 0.0329149

Performing truth solve at parameter:
mu[0] = 2
mu[1] = 0.222222

Enriching the RB space

POD Eigenvalues:
eigenvalue 0 = 0.00014396
eigenvalue 1 = 9.27018e-05
eigenvalue 2 = 2.86671e-05
...
last eigenvalue = -2.10316e-20

Updating RB matrices
Updating RB residual terms
Updating RB initial conditions

---- Basis dimension: 19 ----
Performing RB solves on training set
Maximum (relative) error bound is 0.0256185

Performing truth solve at parameter:
mu[0] = 0.444444
mu[1] = 0

Enriching the RB space

POD Eigenvalues:
eigenvalue 0 = 0.00166579
eigenvalue 1 = 0.000198068
eigenvalue 2 = 0.000175266
...
last eigenvalue = -2.119e-19

Updating RB matrices
Updating RB residual terms
Updating RB initial conditions

---- Basis dimension: 20 ----
Performing RB solves on training set
Maximum (relative) error bound is 0.0234508

Maximum number of basis functions reached: Nmax = 20
Writing out the basis functions...

-------------------------------------------------------------------
| Time:           Sat Apr  7 16:01:42 2012                         |
| OS:             Linux                                            |
| HostName:       lkirk-home                                       |
| OS Release:     3.0.0-17-generic                                 |
| OS Version:     #30-Ubuntu SMP Thu Mar 8 20:45:39 UTC 2012       |
| Machine:        x86_64                                           |
| Username:       benkirk                                          |
| Configuration:  ./configure run on Sat Apr  7 15:49:27 CDT 2012  |
-------------------------------------------------------------------
 -------------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=5.99756, Active time=5.71066                                                      |
 -------------------------------------------------------------------------------------------------------------------
| Event                                 nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                                 w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|-------------------------------------------------------------------------------------------------------------------|
|                                                                                                                   |
|                                                                                                                   |
| DofMap                                                                                                            |
|   add_neighbors_to_send_list()        1         0.0003      0.000341    0.0003      0.000341    0.01     0.01     |
|   build_constraint_matrix()           6875      0.0056      0.000001    0.0056      0.000001    0.10     0.10     |
|   cnstrn_elem_mat_vec()               6875      0.0160      0.000002    0.0160      0.000002    0.28     0.28     |
|   compute_sparsity()                  1         0.0029      0.002868    0.0036      0.003615    0.05     0.06     |
|   create_dof_constraints()            1         0.0026      0.002613    0.0034      0.003443    0.05     0.06     |
|   distribute_dofs()                   1         0.0005      0.000546    0.0021      0.002131    0.01     0.04     |
|   dof_indices()                       15000     0.0111      0.000001    0.0111      0.000001    0.19     0.19     |
|   prepare_send_list()                 1         0.0000      0.000002    0.0000      0.000002    0.00     0.00     |
|   reinit()                            1         0.0016      0.001582    0.0016      0.001582    0.03     0.03     |
|                                                                                                                   |
| EquationSystems                                                                                                   |
|   write()                             1         0.0101      0.010121    0.0104      0.010425    0.18     0.18     |
|                                                                                                                   |
| FE                                                                                                                |
|   compute_affine_map()                7975      0.0059      0.000001    0.0059      0.000001    0.10     0.10     |
|   compute_face_map()                  1100      0.0035      0.000003    0.0072      0.000007    0.06     0.13     |
|   compute_shape_functions()           7975      0.0041      0.000001    0.0041      0.000001    0.07     0.07     |
|   init_face_shape_functions()         11        0.0000      0.000003    0.0000      0.000003    0.00     0.00     |
|   init_shape_functions()              1111      0.0048      0.000004    0.0048      0.000004    0.08     0.08     |
|   inverse_map()                       2200      0.0033      0.000002    0.0033      0.000002    0.06     0.06     |
|                                                                                                                   |
| Mesh                                                                                                              |
|   find_neighbors()                    1         0.0022      0.002228    0.0022      0.002228    0.04     0.04     |
|   renumber_nodes_and_elem()           2         0.0001      0.000071    0.0001      0.000071    0.00     0.00     |
|                                                                                                                   |
| MeshCommunication                                                                                                 |
|   assign_global_indices()             1         0.0161      0.016050    0.0161      0.016062    0.28     0.28     |
|                                                                                                                   |
| MeshTools::Generation                                                                                             |
|   build_cube()                        1         0.0011      0.001140    0.0011      0.001140    0.02     0.02     |
|                                                                                                                   |
| Parallel                                                                                                          |
|   allgather()                         5         0.0000      0.000000    0.0000      0.000000    0.00     0.00     |
|   receive()                           84        0.0003      0.000004    0.0003      0.000004    0.01     0.01     |
|   send()                              84        0.0007      0.000009    0.0007      0.000009    0.01     0.01     |
|   send_receive()                      4         0.0000      0.000002    0.0000      0.000002    0.00     0.00     |
|                                                                                                                   |
| Partitioner                                                                                                       |
|   single_partition()                  1         0.0001      0.000113    0.0001      0.000113    0.00     0.00     |
|                                                                                                                   |
| PetscLinearSolver                                                                                                 |
|   solve()                             2085      1.3768      0.000660    1.3768      0.000660    24.11    24.11    |
|                                                                                                                   |
| RBConstruction                                                                                                    |
|   add_scaled_matrix_and_vector()      11        0.0524      0.004765    0.1137      0.010336    0.92     1.99     |
|   clear()                             3         0.0005      0.000172    0.0005      0.000172    0.01     0.01     |
|   compute_Fq_representor_norms()      1         0.0003      0.000323    0.0034      0.003366    0.01     0.06     |
|   compute_max_error_bound()           21        0.0174      0.000828    1.2890      0.061379    0.30     22.57    |
|   compute_output_dual_norms()         1         0.0105      0.010475    0.0125      0.012476    0.18     0.22     |
|   train_reduced_basis()               1         0.0033      0.003311    5.5917      5.591744    0.06     97.92    |
|   update_RB_system_matrices()         20        0.0545      0.002726    0.0545      0.002726    0.95     0.95     |
|   update_residual_terms()             20        0.1040      0.005199    0.2884      0.014421    1.82     5.05     |
|                                                                                                                   |
| RBEvaluation                                                                                                      |
|   clear()                             2         0.0002      0.000077    0.0003      0.000157    0.00     0.01     |
|   clear_riesz_representors()          4         0.0002      0.000040    0.0002      0.000040    0.00     0.00     |
|   resize_data_structures()            1         0.0001      0.000081    0.0001      0.000082    0.00     0.00     |
|   write_offline_data_to_files()       1         0.0060      0.006048    0.0060      0.006048    0.11     0.11     |
|                                                                                                                   |
| TransientRBConstruction                                                                                           |
|   enrich_RB_space()                   20        0.2968      0.014841    0.2968      0.014841    5.20     5.20     |
|   mass_matrix_scaled_matvec()         2000      0.0994      0.000050    0.0994      0.000050    1.74     1.74     |
|   set_error_temporal_data()           2020      0.2334      0.000116    0.2334      0.000116    4.09     4.09     |
|   truth_assembly()                    2000      1.8915      0.000946    1.9913      0.000996    33.12    34.87    |
|   truth_solve()                       20        0.1424      0.007118    3.5000      0.174998    2.49     61.29    |
|   update_RB_initial_condition_all_N() 20        0.0040      0.000200    0.0040      0.000200    0.07     0.07     |
|   update_RB_system_matrices()         20        0.0176      0.000879    0.0721      0.003606    0.31     1.26     |
|   update_residual_terms()             20        0.0663      0.003313    0.4105      0.020527    1.16     7.19     |
|                                                                                                                   |
| TransientRBEvaluation                                                                                             |
|   cache_online_residual_terms()       2100      0.0083      0.000004    0.0083      0.000004    0.15     0.15     |
|   compute_residual_dual_norm()        210000    0.3416      0.000002    0.3416      0.000002    5.98     5.98     |
|   rb_solve()                          2100      0.8866      0.000422    1.2712      0.000605    15.52    22.26    |
|   resize_data_structures()            1         0.0001      0.000055    0.0001      0.000137    0.00     0.00     |
|   write_offline_data_to_files()       1         0.0035      0.003453    0.0095      0.009501    0.06     0.17     |
 -------------------------------------------------------------------------------------------------------------------
| Totals:                               271805    5.7107                                          100.00            |
 -------------------------------------------------------------------------------------------------------------------

*** Warning, This code is untested, experimental, or likely to see future API changes: src/reduced_basis/rb_parametrized_object.C, line 31, compiled Apr  7 2012 at 15:51:42 ***
mu[0] = 1
mu[1] = 0.79

Error bound (absolute) at the final time is 0.00525183

Reading in the basis functions...
Finished reading in the basis functions...

-------------------------------------------------------------------
| Time:           Sat Apr  7 16:01:42 2012                         |
| OS:             Linux                                            |
| HostName:       lkirk-home                                       |
| OS Release:     3.0.0-17-generic                                 |
| OS Version:     #30-Ubuntu SMP Thu Mar 8 20:45:39 UTC 2012       |
| Machine:        x86_64                                           |
| Username:       benkirk                                          |
| Configuration:  ./configure run on Sat Apr  7 15:49:27 CDT 2012  |
-------------------------------------------------------------------
 --------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=0.266002, Active time=0.056725                                               |
 --------------------------------------------------------------------------------------------------------------
| Event                            nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                            w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|--------------------------------------------------------------------------------------------------------------|
|                                                                                                              |
|                                                                                                              |
| DofMap                                                                                                       |
|   add_neighbors_to_send_list()   1         0.0003      0.000340    0.0003      0.000340    0.60     0.60     |
|   compute_sparsity()             1         0.0020      0.002014    0.0026      0.002615    3.55     4.61     |
|   create_dof_constraints()       1         0.0003      0.000323    0.0003      0.000323    0.57     0.57     |
|   distribute_dofs()              1         0.0005      0.000540    0.0021      0.002105    0.95     3.71     |
|   dof_indices()                  1875      0.0012      0.000001    0.0012      0.000001    2.17     2.17     |
|   prepare_send_list()            1         0.0000      0.000001    0.0000      0.000001    0.00     0.00     |
|   reinit()                       1         0.0016      0.001564    0.0016      0.001564    2.76     2.76     |
|                                                                                                              |
| EquationSystems                                                                                              |
|   build_solution_vector()        2         0.0016      0.000779    0.0025      0.001269    2.75     4.47     |
|   read()                         1         0.0113      0.011334    0.0305      0.030511    19.98    53.79    |
|   update()                       1         0.0001      0.000051    0.0001      0.000051    0.09     0.09     |
|                                                                                                              |
| ExodusII_IO                                                                                                  |
|   write_nodal_data()             2         0.0029      0.001473    0.0029      0.001473    5.19     5.19     |
|                                                                                                              |
| Mesh                                                                                                         |
|   find_neighbors()               1         0.0021      0.002102    0.0021      0.002102    3.71     3.71     |
|   renumber_nodes_and_elem()      2         0.0001      0.000066    0.0001      0.000066    0.23     0.23     |
|                                                                                                              |
| MeshCommunication                                                                                            |
|   assign_global_indices()        1         0.0137      0.013718    0.0137      0.013726    24.18    24.20    |
|                                                                                                              |
| MeshOutput                                                                                                   |
|   write_equation_systems()       2         0.0000      0.000022    0.0055      0.002766    0.08     9.75     |
|                                                                                                              |
| MeshTools::Generation                                                                                        |
|   build_cube()                   1         0.0012      0.001152    0.0012      0.001152    2.03     2.03     |
|                                                                                                              |
| Parallel                                                                                                     |
|   allgather()                    47        0.0001      0.000001    0.0001      0.000001    0.11     0.11     |
|   send_receive()                 4         0.0000      0.000001    0.0000      0.000001    0.01     0.01     |
|                                                                                                              |
| Partitioner                                                                                                  |
|   single_partition()             1         0.0001      0.000134    0.0001      0.000134    0.24     0.24     |
|                                                                                                              |
| RBConstruction                                                                                               |
|   clear()                        4         0.0002      0.000056    0.0002      0.000056    0.40     0.40     |
|   load_basis_function()          1         0.0000      0.000024    0.0000      0.000024    0.04     0.04     |
|                                                                                                              |
| RBEvaluation                                                                                                 |
|   clear()                        2         0.0000      0.000021    0.0000      0.000023    0.08     0.08     |
|   clear_riesz_representors()     4         0.0000      0.000002    0.0000      0.000002    0.01     0.01     |
|   read_offline_data_from_files() 1         0.0089      0.008852    0.0092      0.009178    15.61    16.18    |
|   resize_data_structures()       1         0.0002      0.000209    0.0002      0.000214    0.37     0.38     |
|                                                                                                              |
| TransientRBConstruction                                                                                      |
|   load_rb_solution()             1         0.0001      0.000120    0.0001      0.000120    0.21     0.21     |
|                                                                                                              |
| TransientRBEvaluation                                                                                        |
|   cache_online_residual_terms()  1         0.0000      0.000036    0.0000      0.000036    0.06     0.06     |
|   compute_residual_dual_norm()   100       0.0008      0.000008    0.0008      0.000008    1.34     1.34     |
|   rb_solve()                     1         0.0019      0.001862    0.0027      0.002689    3.28     4.74     |
|   read_offline_data_from_files() 1         0.0052      0.005222    0.0144      0.014400    9.21     25.39    |
|   resize_data_structures()       1         0.0001      0.000111    0.0003      0.000325    0.20     0.57     |
 --------------------------------------------------------------------------------------------------------------
| Totals:                          2064      0.0567                                          100.00            |
 --------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running  ./reduced_basis_ex3-opt
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
