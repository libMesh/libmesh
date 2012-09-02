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
  

<br><br></div>

<div class ="fragment">
<pre>
        /* rbOOmit is distributed in the hope that it will be useful, */
        /* but WITHOUT ANY WARRANTY; without even the implied warranty of */
        /* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU */
        /* Lesser General Public License for more details. */
          
        /* You should have received a copy of the GNU Lesser General Public */
        /* License along with this library; if not, write to the Free Software */
        /* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA */
        
</pre>
</div>
<div class = "comment">
<h1>Reduced Basis Example 3 - Transient Reduced Basis Problem</h1>


<br><br>In this example problem we use the Certified Reduced Basis method
to solve a transient convection-diffusion problem on the unit square.
The PDE is similar to reduced_basis_ex1, except there is a time-derivative
in this case.


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
Build a new RBEvaluation object which will be used to perform
Reduced Basis calculations. This is required in both the
"Offline" and "Online" stages.
</div>

<div class ="fragment">
<pre>
          SimpleRBEvaluation rb_eval;
        
</pre>
</div>
<div class = "comment">
Finally, we need to give the RBConstruction object a pointer to
our RBEvaluation object
</div>

<div class ="fragment">
<pre>
          rb_con.set_rb_evaluation(rb_eval);
        
          if(!online_mode) // Perform the Offline stage of the RB method
          {
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
            rb_con.get_rb_evaluation().write_offline_data_to_files();
            
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
Write out the basis functions
</div>

<div class ="fragment">
<pre>
              rb_con.get_rb_evaluation().write_out_basis_functions(rb_con);
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
Read in online_N and initialize online parameters
</div>

<div class ="fragment">
<pre>
            unsigned int online_N = infile("online_N",1);
            Real online_x_vel = infile("online_x_vel", 0.);
            Real online_y_vel = infile("online_y_vel", 0.);
            RBParameters online_mu;
            online_mu.set_value("x_vel", online_x_vel);
            online_mu.set_value("y_vel", online_y_vel);
            rb_eval.set_parameters(online_mu);
            rb_eval.print_parameters();
        
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
              rb_con.pull_temporal_discretization_data( rb_eval );
              rb_con.set_time_step(rb_con.get_n_time_steps());
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
    
  <I><FONT COLOR="#B22222">/* rbOOmit is distributed in the hope that it will be useful, */</FONT></I>
  <I><FONT COLOR="#B22222">/* but WITHOUT ANY WARRANTY; without even the implied warranty of */</FONT></I>
  <I><FONT COLOR="#B22222">/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU */</FONT></I>
  <I><FONT COLOR="#B22222">/* Lesser General Public License for more details. */</FONT></I>
    
  <I><FONT COLOR="#B22222">/* You should have received a copy of the GNU Lesser General Public */</FONT></I>
  <I><FONT COLOR="#B22222">/* License along with this library; if not, write to the Free Software */</FONT></I>
  <I><FONT COLOR="#B22222">/* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA */</FONT></I>
  
  
  
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
    
    SimpleRBConstruction &amp; rb_con =
      equation_systems.add_system&lt;SimpleRBConstruction&gt; (<B><FONT COLOR="#BC8F8F">&quot;RBConvectionDiffusion&quot;</FONT></B>);
  
  
    equation_systems.init ();
  
    equation_systems.print_info();
    mesh.print_info();
  
    SimpleRBEvaluation rb_eval;
  
    rb_con.set_rb_evaluation(rb_eval);
  
    <B><FONT COLOR="#A020F0">if</FONT></B>(!online_mode) <I><FONT COLOR="#B22222">// Perform the Offline stage of the RB method
</FONT></I>    {
      rb_con.process_parameters_file(parameters_filename);
  
      rb_con.print_info();
  
      rb_con.initialize_rb_construction();
  
      rb_con.train_reduced_basis();
      
      rb_con.get_rb_evaluation().write_offline_data_to_files();
      
      <B><FONT COLOR="#A020F0">if</FONT></B>(store_basis_functions)
      {
        rb_con.get_rb_evaluation().write_out_basis_functions(rb_con);
      }
    }
    <B><FONT COLOR="#A020F0">else</FONT></B> <I><FONT COLOR="#B22222">// Perform the Online stage of the RB method
</FONT></I>    {
      rb_eval.read_offline_data_from_files();
      
      <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> online_N = infile(<B><FONT COLOR="#BC8F8F">&quot;online_N&quot;</FONT></B>,1);
      Real online_x_vel = infile(<B><FONT COLOR="#BC8F8F">&quot;online_x_vel&quot;</FONT></B>, 0.);
      Real online_y_vel = infile(<B><FONT COLOR="#BC8F8F">&quot;online_y_vel&quot;</FONT></B>, 0.);
      RBParameters online_mu;
      online_mu.set_value(<B><FONT COLOR="#BC8F8F">&quot;x_vel&quot;</FONT></B>, online_x_vel);
      online_mu.set_value(<B><FONT COLOR="#BC8F8F">&quot;y_vel&quot;</FONT></B>, online_y_vel);
      rb_eval.set_parameters(online_mu);
      rb_eval.print_parameters();
  
      Real error_bound_final_time = rb_eval.rb_solve(online_N);
      
      <B><FONT COLOR="#5F9EA0">libMesh</FONT></B>::out &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Error bound (absolute) at the final time is &quot;</FONT></B>
                   &lt;&lt; error_bound_final_time &lt;&lt; std::endl &lt;&lt; std::endl;
  
      <B><FONT COLOR="#A020F0">if</FONT></B>(store_basis_functions)
      {
        rb_eval.read_in_basis_functions(rb_con);
        
        rb_con.pull_temporal_discretization_data( rb_eval );
        rb_con.set_time_step(rb_con.get_n_time_steps());
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
Linking reduced_basis_ex3-opt...
***************************************************************
* Running  mpirun -np 6 ./reduced_basis_ex3-opt
***************************************************************
 
*** Warning, This code is untested, experimental, or likely to see future API changes: src/reduced_basis/rb_parametrized.C, line 40, compiled Aug 24 2012 at 15:15:42 ***
 EquationSystems
  n_systems()=1
   System #0, "RBConvectionDiffusion"
    Type "TransientRBConstruction"
    Variables="u" 
    Finite Element Types="LAGRANGE", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" 
    Approximation Orders="FIRST", "THIRD" 
    n_dofs()=676
    n_local_dofs()=129
    n_constrained_dofs()=102
    n_local_constrained_dofs()=21
    n_vectors()=3
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 7.94574
      Average Off-Processor Bandwidth <= 0.550388
      Maximum  On-Processor Bandwidth <= 9
      Maximum Off-Processor Bandwidth <= 5
    DofMap Constraints
      Number of DoF Constraints = 100
      Average DoF Constraint Length= 0
      Number of Node Constraints = 0

 Mesh Information:
  mesh_dimension()=2
  spatial_dimension()=3
  n_nodes()=676
    n_local_nodes()=129
  n_elem()=625
    n_local_elem()=106
    n_active_elem()=625
  n_subdomains()=1
  n_partitions()=6
  n_processors()=6
  n_threads()=1
  processor_id()=0

Initializing training parameters with deterministic training set...
Parameter x_vel: log scaling = 0
Parameter y_vel: log scaling = 0


RBConstruction parameters:
system name: RBConvectionDiffusion
constrained_problem: 0
Nmax: 20
Aq operators attached: 3
Fq functions attached: 1
n_outputs: 4
output 0, Q_l = 1
output 1, Q_l = 1
output 2, Q_l = 1
output 3, Q_l = 1
Number of parameters: 2
Parameter x_vel: Min = -2, Max = 2, value = 1
Parameter y_vel: Min = -2, Max = 2, value = 1
n_training_samples: 100
single-matrix mode? 0
reuse preconditioner? 1
use a relative error bound in greedy? 1
write out data during basis training? 0
quiet mode? 1


TransientRBConstruction parameters:
Q_m: 1
Number of time-steps: 100
dt: 0.01
euler_theta (time discretization parameter): 1
delta_N (number of basis functions to add each POD-Greedy step): 1
Using zero initial condition

Compute output dual norms
output_dual_innerprods[0][0] = 33.6538
output_dual_innerprods[1][0] = 33.6538
output_dual_innerprods[2][0] = 33.6538
output_dual_innerprods[3][0] = 33.6538

---- Performing Greedy basis enrichment ----

---- Basis dimension: 0 ----
Performing RB solves on training set
Maximum (relative) error bound is inf

Performing truth solve at parameter:
x_vel: -2
y_vel: -2

Enriching the RB space

POD Eigenvalues:
eigenvalue 0 = 2.19193
eigenvalue 1 = 0.0116709
eigenvalue 2 = 0.00106473
...
last eigenvalue = -2.84727e-16

Updating RB matrices
Updating RB residual terms
Updating RB initial conditions

---- Basis dimension: 1 ----
Performing RB solves on training set
Maximum (relative) error bound is 6.81811

Performing truth solve at parameter:
x_vel: 2
y_vel: 2

Enriching the RB space

POD Eigenvalues:
eigenvalue 0 = 2.14617
eigenvalue 1 = 0.0106903
eigenvalue 2 = 0.000908038
...
last eigenvalue = -3.37378e-16

Updating RB matrices
Updating RB residual terms
Updating RB initial conditions

---- Basis dimension: 2 ----
Performing RB solves on training set
Maximum (relative) error bound is 3.15659

Performing truth solve at parameter:
x_vel: 2
y_vel: -2

Enriching the RB space

POD Eigenvalues:
eigenvalue 0 = 1.6024
eigenvalue 1 = 0.00595727
eigenvalue 2 = 0.000820407
...
last eigenvalue = -1.61648e-16

Updating RB matrices
Updating RB residual terms
Updating RB initial conditions

---- Basis dimension: 3 ----
Performing RB solves on training set
Maximum (relative) error bound is 3.95604

Performing truth solve at parameter:
x_vel: -2
y_vel: 2

Enriching the RB space

POD Eigenvalues:
eigenvalue 0 = 1.55558
eigenvalue 1 = 0.0058283
eigenvalue 2 = 0.00079201
...
last eigenvalue = -6.37081e-17

Updating RB matrices
Updating RB residual terms
Updating RB initial conditions

---- Basis dimension: 4 ----
Performing RB solves on training set
Maximum (relative) error bound is 1.3165

Performing truth solve at parameter:
x_vel: -0.222222
y_vel: 0.222222

Enriching the RB space

POD Eigenvalues:
eigenvalue 0 = 1.02449
eigenvalue 1 = 0.0148056
eigenvalue 2 = 0.000593638
...
last eigenvalue = -2.17328e-16

Updating RB matrices
Updating RB residual terms
Updating RB initial conditions

---- Basis dimension: 5 ----
Performing RB solves on training set
Maximum (relative) error bound is 1.05097

Performing truth solve at parameter:
x_vel: -2
y_vel: -0.222222

Enriching the RB space

POD Eigenvalues:
eigenvalue 0 = 0.365693
eigenvalue 1 = 0.00367148
eigenvalue 2 = 0.000675947
...
last eigenvalue = -7.51469e-18

Updating RB matrices
Updating RB residual terms
Updating RB initial conditions

---- Basis dimension: 6 ----
Performing RB solves on training set
Maximum (relative) error bound is 1.0171

Performing truth solve at parameter:
x_vel: -0.222222
y_vel: -2

Enriching the RB space

POD Eigenvalues:
eigenvalue 0 = 0.285715
eigenvalue 1 = 0.00305589
eigenvalue 2 = 0.000804404
...
last eigenvalue = -1.44343e-18

Updating RB matrices
Updating RB residual terms
Updating RB initial conditions

---- Basis dimension: 7 ----
Performing RB solves on training set
Maximum (relative) error bound is 1.14302

Performing truth solve at parameter:
x_vel: -0.222222
y_vel: 2

Enriching the RB space

POD Eigenvalues:
eigenvalue 0 = 0.256366
eigenvalue 1 = 0.00345875
eigenvalue 2 = 0.000587394
...
last eigenvalue = -1.68438e-17

Updating RB matrices
Updating RB residual terms
Updating RB initial conditions

---- Basis dimension: 8 ----
Performing RB solves on training set
Maximum (relative) error bound is 1.21494

Performing truth solve at parameter:
x_vel: 2
y_vel: -0.222222

Enriching the RB space

POD Eigenvalues:
eigenvalue 0 = 0.260659
eigenvalue 1 = 0.00277678
eigenvalue 2 = 0.000773864
...
last eigenvalue = -1.08439e-18

Updating RB matrices
Updating RB residual terms
Updating RB initial conditions

---- Basis dimension: 9 ----
Performing RB solves on training set
Maximum (relative) error bound is 0.43786

Performing truth solve at parameter:
x_vel: 1.11111
y_vel: 1.11111

Enriching the RB space

POD Eigenvalues:
eigenvalue 0 = 0.0933151
eigenvalue 1 = 0.00402501
eigenvalue 2 = 0.00138016
...
last eigenvalue = -3.92937e-18

Updating RB matrices
Updating RB residual terms
Updating RB initial conditions

---- Basis dimension: 10 ----
Performing RB solves on training set
Maximum (relative) error bound is 0.433787

Performing truth solve at parameter:
x_vel: 0.666667
y_vel: -1.11111

Enriching the RB space

POD Eigenvalues:
eigenvalue 0 = 0.134913
eigenvalue 1 = 0.00436435
eigenvalue 2 = 0.00137045
...
last eigenvalue = -8.85445e-18

Updating RB matrices
Updating RB residual terms
Updating RB initial conditions

---- Basis dimension: 11 ----
Performing RB solves on training set
Maximum (relative) error bound is 0.319326

Performing truth solve at parameter:
x_vel: -1.11111
y_vel: -1.11111

Enriching the RB space

POD Eigenvalues:
eigenvalue 0 = 0.0629291
eigenvalue 1 = 0.00315345
eigenvalue 2 = 0.00120446
...
last eigenvalue = -1.60013e-18

Updating RB matrices
Updating RB residual terms
Updating RB initial conditions

---- Basis dimension: 12 ----
Performing RB solves on training set
Maximum (relative) error bound is 0.287715

Performing truth solve at parameter:
x_vel: -1.55556
y_vel: 1.11111

Enriching the RB space

POD Eigenvalues:
eigenvalue 0 = 0.0278364
eigenvalue 1 = 0.00277066
eigenvalue 2 = 0.000760673
...
last eigenvalue = -1.1209e-19

Updating RB matrices
Updating RB residual terms
Updating RB initial conditions

---- Basis dimension: 13 ----
Performing RB solves on training set
Maximum (relative) error bound is 0.228017

Performing truth solve at parameter:
x_vel: 1.11111
y_vel: -2

Enriching the RB space

POD Eigenvalues:
eigenvalue 0 = 0.00863589
eigenvalue 1 = 0.00176423
eigenvalue 2 = 0.000525486
...
last eigenvalue = -1.02632e-18

Updating RB matrices
Updating RB residual terms
Updating RB initial conditions

---- Basis dimension: 14 ----
Performing RB solves on training set
Maximum (relative) error bound is 0.227118

Performing truth solve at parameter:
x_vel: 2
y_vel: 1.11111

Enriching the RB space

POD Eigenvalues:
eigenvalue 0 = 0.0117352
eigenvalue 1 = 0.00206357
eigenvalue 2 = 0.000546886
...
last eigenvalue = -6.52666e-19

Updating RB matrices
Updating RB residual terms
Updating RB initial conditions

---- Basis dimension: 15 ----
Performing RB solves on training set
Maximum (relative) error bound is 0.231987

Performing truth solve at parameter:
x_vel: 0.666667
y_vel: 2

Enriching the RB space

POD Eigenvalues:
eigenvalue 0 = 0.00898933
eigenvalue 1 = 0.0022199
eigenvalue 2 = 0.000417343
...
last eigenvalue = -9.57807e-19

Updating RB matrices
Updating RB residual terms
Updating RB initial conditions

---- Basis dimension: 16 ----
Performing RB solves on training set
Maximum (relative) error bound is 0.199807

Performing truth solve at parameter:
x_vel: -2
y_vel: 0.666667

Enriching the RB space

POD Eigenvalues:
eigenvalue 0 = 0.0065467
eigenvalue 1 = 0.00176326
eigenvalue 2 = 0.000231436
...
last eigenvalue = -5.55006e-19

Updating RB matrices
Updating RB residual terms
Updating RB initial conditions

---- Basis dimension: 17 ----
Performing RB solves on training set
Maximum (relative) error bound is 0.193551

Performing truth solve at parameter:
x_vel: -1.11111
y_vel: 2

Enriching the RB space

POD Eigenvalues:
eigenvalue 0 = 0.00390004
eigenvalue 1 = 0.00190681
eigenvalue 2 = 0.000456244
...
last eigenvalue = -6.09332e-20

Updating RB matrices
Updating RB residual terms
Updating RB initial conditions

---- Basis dimension: 18 ----
Performing RB solves on training set
Maximum (relative) error bound is 0.181413

Performing truth solve at parameter:
x_vel: -1.11111
y_vel: -2

Enriching the RB space

POD Eigenvalues:
eigenvalue 0 = 0.00391101
eigenvalue 1 = 0.00171004
eigenvalue 2 = 0.000346167
...
last eigenvalue = -3.65918e-19

Updating RB matrices
Updating RB residual terms
Updating RB initial conditions

---- Basis dimension: 19 ----
Performing RB solves on training set
Maximum (relative) error bound is 0.174248

Performing truth solve at parameter:
x_vel: 1.11111
y_vel: -0.222222

Enriching the RB space

POD Eigenvalues:
eigenvalue 0 = 0.0287899
eigenvalue 1 = 0.00198823
eigenvalue 2 = 0.000531507
...
last eigenvalue = -3.69478e-19

Updating RB matrices
Updating RB residual terms
Updating RB initial conditions

---- Basis dimension: 20 ----
Performing RB solves on training set
Maximum (relative) error bound is 0.145574

Maximum number of basis functions reached: Nmax = 20
************************************************************************************************************************
***             WIDEN YOUR WINDOW TO 120 CHARACTERS.  Use 'enscript -r -fCourier9' to print this document            ***
************************************************************************************************************************

---------------------------------------------- PETSc Performance Summary: ----------------------------------------------

./reduced_basis_ex3-opt on a intel-11. named daedalus with 6 processors, by roystgnr Fri Aug 24 15:22:40 2012
Using Petsc Release Version 3.1.0, Patch 5, Mon Sep 27 11:51:54 CDT 2010

                         Max       Max/Min        Avg      Total 
Time (sec):           1.826e+01      1.00006   1.826e+01
Objects:              6.502e+03      1.00000   6.502e+03
Flops:                3.553e+08      1.33701   3.036e+08  1.822e+09
Flops/sec:            1.946e+07      1.33708   1.663e+07  9.976e+07
MPI Messages:         2.469e+05      2.50000   1.488e+05  8.931e+05
MPI Message Lengths:  2.023e+07      2.06703   8.978e+01  8.018e+07
MPI Reductions:       3.470e+05      1.00000

Flop counting convention: 1 flop = 1 real number operation of type (multiply/divide/add/subtract)
                            e.g., VecAXPY() for real vectors of length N --> 2N flops
                            and VecAXPY() for complex vectors of length N --> 8N flops

Summary of Stages:   ----- Time ------  ----- Flops -----  --- Messages ---  -- Message Lengths --  -- Reductions --
                        Avg     %Total     Avg     %Total   counts   %Total     Avg         %Total   counts   %Total 
 0:      Main Stage: 1.8258e+01 100.0%  1.8215e+09 100.0%  8.931e+05 100.0%  8.978e+01      100.0%  3.407e+05  98.2% 

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

VecDot            138405 1.0 5.6478e+00 1.1 3.56e+07 1.3 0.0e+00 0.0e+00 1.4e+05 30 10  0  0 40  30 10  0  0 41    33
VecMDot            23066 1.0 1.1820e+00 1.9 4.06e+07 1.3 0.0e+00 0.0e+00 2.3e+04  5 12  0  0  7   5 12  0  0  7   180
VecNorm            27317 1.0 1.3611e+00 1.2 7.05e+06 1.3 0.0e+00 0.0e+00 2.7e+04  7  2  0  0  8   7  2  0  0  8    27
VecScale           33312 1.0 1.5764e-02 1.3 3.27e+06 1.3 0.0e+00 0.0e+00 0.0e+00  0  1  0  0  0   0  1  0  0  0  1085
VecCopy            14627 1.0 5.2910e-03 1.3 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecSet             61734 1.0 2.2749e-02 1.3 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecAXPY            37797 1.0 1.7025e-02 1.3 9.64e+06 1.3 0.0e+00 0.0e+00 0.0e+00  0  3  0  0  0   0  3  0  0  0  2968
VecMAXPY           25232 1.0 2.4590e-02 1.3 4.67e+07 1.3 0.0e+00 0.0e+00 0.0e+00  0 13  0  0  0   0 13  0  0  0  9958
VecAssemblyBegin   14344 1.0 1.8122e+00 1.0 0.00e+00 0.0 9.0e+01 1.3e+02 4.3e+04 10  0  0  0 12  10  0  0  0 13     0
VecAssemblyEnd     14344 1.0 8.4329e-03 1.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecScatterBegin    49351 1.0 1.2725e-01 1.7 0.00e+00 0.0 8.9e+05 9.0e+01 0.0e+00  1  0100100  0   1  0100100  0     0
VecScatterEnd      49351 1.0 2.1970e+00 2.7 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  7  0  0  0  0   7  0  0  0  0     0
VecNormalize       25232 1.0 1.2747e+00 1.2 9.76e+06 1.3 0.0e+00 0.0e+00 2.5e+04  6  3  0  0  7   6  3  0  0  7    40
MatMult            25232 1.0 1.2836e+00 2.3 5.21e+07 1.3 4.5e+05 8.6e+01 0.0e+00  5 15 51 49  0   5 15 51 49  0   214
MatMultAdd         19994 1.0 1.0164e+00 2.6 4.38e+07 1.3 3.6e+05 8.6e+01 0.0e+00  3 13 40 39  0   3 13 40 39  0   227
MatSolve           25232 1.0 1.0098e-01 1.5 1.16e+08 1.6 0.0e+00 0.0e+00 0.0e+00  0 31  0  0  0   0 31  0  0  0  5588
MatLUFactorNum        42 1.0 1.9886e-03 1.4 7.53e+05 2.2 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0  1708
MatILUFactorSym       42 1.0 4.8993e-03 1.7 0.00e+00 0.0 0.0e+00 0.0e+00 4.2e+01  0  0  0  0  0   0  0  0  0  0     0
MatAssemblyBegin   36215 1.0 3.7851e+00 1.1 0.00e+00 0.0 1.6e+02 4.3e+02 7.2e+04 20  0  0  0 21  20  0  0  0 21     0
MatAssemblyEnd     36215 1.0 1.7949e+00 1.1 0.00e+00 0.0 2.5e+02 2.4e+01 3.6e+04  9  0  0  0 10   9  0  0  0 11     0
MatGetRow        1037418 1.3 1.1867e-01 1.3 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  1  0  0  0  0   1  0  0  0  0     0
MatGetRowIJ           42 1.0 1.4544e-05 1.8 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatGetOrdering        42 1.0 3.5453e-04 1.2 0.00e+00 0.0 0.0e+00 0.0e+00 8.4e+01  0  0  0  0  0   0  0  0  0  0     0
MatZeroEntries      2062 1.0 2.6507e-03 1.5 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPGMRESOrthog     23066 1.0 1.2124e+00 1.9 8.14e+07 1.3 0.0e+00 0.0e+00 2.3e+04  6 23  0  0  7   6 23  0  0  7   352
KSPSetup            2127 1.0 3.2759e-04 1.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPSolve            2085 1.0 3.3970e+00 1.0 2.67e+08 1.4 4.5e+05 8.6e+01 5.1e+04 19 75 51 49 15  19 75 51 49 15   400
PCSetUp               84 1.0 9.1434e-03 1.5 7.53e+05 2.2 0.0e+00 0.0e+00 1.3e+02  0  0  0  0  0   0  0  0  0  0   372
PCSetUpOnBlocks     2085 1.0 8.9903e-03 1.4 7.53e+05 2.2 0.0e+00 0.0e+00 1.3e+02  0  0  0  0  0   0  0  0  0  0   378
PCApply            25232 1.0 2.5553e-01 1.2 1.16e+08 1.6 0.0e+00 0.0e+00 0.0e+00  1 31  0  0  0   1 31  0  0  0  2208
------------------------------------------------------------------------------------------------------------------------

Memory usage is given in bytes:

Object Type          Creations   Destructions     Memory  Descendants' Mem.
Reports information only for process 0.

--- Event Stage 0: Main Stage

                 Vec  6272           6272     14880136     0
         Vec Scatter    11             11         9548     0
           Index Set   148            148       143032     0
   IS L to G Mapping     3              3         3084     0
              Matrix    64             64      1454824     0
       Krylov Solver     2              2        18880     0
      Preconditioner     2              2         1408     0
========================================================================================================================
Average time to get PetscTime(): 0
Average time for MPI_Barrier(): 4.00066e-05
Average time for zero size MPI_Send(): 4.61737e-05
#PETSc Option Table entries:
-ksp_right_pc
-log_summary
-online_mode 0
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
| Time:           Fri Aug 24 15:22:40 2012                         |
| OS:             Linux                                            |
| HostName:       daedalus                                         |
| OS Release:     2.6.32-34-generic                                |
| OS Version:     #76-Ubuntu SMP Tue Aug 30 17:05:01 UTC 2011      |
| Machine:        x86_64                                           |
| Username:       roystgnr                                         |
| Configuration:  ./configure run on Wed Aug 22 12:44:06 CDT 2012  |
-------------------------------------------------------------------
 -------------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=18.4377, Active time=18.225                                                       |
 -------------------------------------------------------------------------------------------------------------------
| Event                                 nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                                 w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|-------------------------------------------------------------------------------------------------------------------|
|                                                                                                                   |
|                                                                                                                   |
| DofMap                                                                                                            |
|   add_neighbors_to_send_list()        1         0.0013      0.001300    0.0014      0.001367    0.01     0.01     |
|   build_constraint_matrix()           1166      0.0008      0.000001    0.0008      0.000001    0.00     0.00     |
|   build_sparsity()                    1         0.0006      0.000648    0.0008      0.000849    0.00     0.00     |
|   cnstrn_elem_mat_vec()               1166      0.0079      0.000007    0.0079      0.000007    0.04     0.04     |
|   create_dof_constraints()            1         0.0008      0.000831    0.0010      0.001024    0.00     0.01     |
|   distribute_dofs()                   1         0.0002      0.000206    0.0013      0.001306    0.00     0.01     |
|   dof_indices()                       3631      0.0012      0.000000    0.0012      0.000000    0.01     0.01     |
|   prepare_send_list()                 1         0.0000      0.000003    0.0000      0.000003    0.00     0.00     |
|   reinit()                            1         0.0004      0.000398    0.0004      0.000398    0.00     0.00     |
|                                                                                                                   |
| FE                                                                                                                |
|   compute_shape_functions()           2772      0.0037      0.000001    0.0037      0.000001    0.02     0.02     |
|   init_shape_functions()              462       0.0014      0.000003    0.0014      0.000003    0.01     0.01     |
|   inverse_map()                       880       0.0008      0.000001    0.0008      0.000001    0.00     0.00     |
|                                                                                                                   |
| FEMap                                                                                                             |
|   compute_affine_map()                2772      0.0017      0.000001    0.0017      0.000001    0.01     0.01     |
|   compute_face_map()                  440       0.0010      0.000002    0.0019      0.000004    0.01     0.01     |
|   init_face_shape_functions()         22        0.0000      0.000002    0.0000      0.000002    0.00     0.00     |
|   init_reference_to_physical_map()    462       0.0008      0.000002    0.0008      0.000002    0.00     0.00     |
|                                                                                                                   |
| Mesh                                                                                                              |
|   find_neighbors()                    1         0.0005      0.000550    0.0057      0.005690    0.00     0.03     |
|   renumber_nodes_and_elem()           2         0.0000      0.000020    0.0000      0.000020    0.00     0.00     |
|                                                                                                                   |
| MeshCommunication                                                                                                 |
|   assign_global_indices()             1         0.0057      0.005650    0.0104      0.010374    0.03     0.06     |
|   compute_hilbert_indices()           2         0.0036      0.001802    0.0036      0.001802    0.02     0.02     |
|   find_global_indices()               2         0.0002      0.000121    0.0077      0.003850    0.00     0.04     |
|   parallel_sort()                     2         0.0015      0.000758    0.0030      0.001519    0.01     0.02     |
|                                                                                                                   |
| MeshTools::Generation                                                                                             |
|   build_cube()                        1         0.0003      0.000260    0.0003      0.000260    0.00     0.00     |
|                                                                                                                   |
| MetisPartitioner                                                                                                  |
|   partition()                         1         0.0028      0.002805    0.0063      0.006311    0.02     0.03     |
|                                                                                                                   |
| Parallel                                                                                                          |
|   allgather()                         13        0.0045      0.000345    0.0045      0.000345    0.02     0.02     |
|   barrier()                           21        0.0007      0.000031    0.0007      0.000031    0.00     0.00     |
|   broadcast()                         43        0.0011      0.000026    0.0011      0.000026    0.01     0.01     |
|   gather()                            81        0.0035      0.000044    0.0035      0.000044    0.02     0.02     |
|   max(scalar)                         2         0.0051      0.002572    0.0051      0.002572    0.03     0.03     |
|   max(vector)                         3         0.0003      0.000107    0.0003      0.000107    0.00     0.00     |
|   maxloc(scalar)                      21        0.1365      0.006499    0.1365      0.006499    0.75     0.75     |
|   min(vector)                         3         0.0006      0.000197    0.0006      0.000197    0.00     0.00     |
|   probe()                             70        0.0011      0.000015    0.0011      0.000015    0.01     0.01     |
|   receive()                           550       0.0005      0.000001    0.0015      0.000003    0.00     0.01     |
|   send()                              150       0.0001      0.000001    0.0001      0.000001    0.00     0.00     |
|   send_receive()                      78        0.0001      0.000002    0.0014      0.000018    0.00     0.01     |
|   sum()                               34        0.0087      0.000257    0.0087      0.000257    0.05     0.05     |
|                                                                                                                   |
| Parallel::Request                                                                                                 |
|   wait()                              550       0.0003      0.000000    0.0003      0.000000    0.00     0.00     |
|                                                                                                                   |
| Partitioner                                                                                                       |
|   set_node_processor_ids()            1         0.0001      0.000081    0.0002      0.000214    0.00     0.00     |
|   set_parent_processor_ids()          1         0.0000      0.000039    0.0000      0.000039    0.00     0.00     |
|                                                                                                                   |
| PetscLinearSolver                                                                                                 |
|   solve()                             2085      4.4130      0.002117    4.4130      0.002117    24.21    24.21    |
|                                                                                                                   |
| RBConstruction                                                                                                    |
|   add_scaled_matrix_and_vector()      11        0.0420      0.003821    0.0618      0.005619    0.23     0.34     |
|   clear()                             3         0.0003      0.000101    0.0003      0.000101    0.00     0.00     |
|   compute_Fq_representor_innerprods() 1         0.0012      0.001187    0.0125      0.012548    0.01     0.07     |
|   compute_max_error_bound()           21        0.0022      0.000104    0.2845      0.013549    0.01     1.56     |
|   compute_output_dual_innerprods()    1         0.0054      0.005356    0.0221      0.022068    0.03     0.12     |
|   train_reduced_basis()               1         0.0246      0.024591    18.1095     18.109468   0.13     99.37    |
|   update_RB_system_matrices()         20        0.2936      0.014680    0.2936      0.014680    1.61     1.61     |
|   update_residual_terms()             20        0.5683      0.028416    1.0778      0.053890    3.12     5.91     |
|                                                                                                                   |
| RBEvaluation                                                                                                      |
|   clear()                             2         0.0001      0.000045    0.0001      0.000045    0.00     0.00     |
|   resize_data_structures()            1         0.0001      0.000076    0.0001      0.000076    0.00     0.00     |
|   write_offline_data_to_files()       1         0.0197      0.019726    0.0197      0.019726    0.11     0.11     |
|                                                                                                                   |
| TransientRBConstruction                                                                                           |
|   enrich_RB_space()                   20        4.4162      0.220809    4.4162      0.220809    24.23    24.23    |
|   mass_matrix_scaled_matvec()         2000      0.4589      0.000229    0.4589      0.000229    2.52     2.52     |
|   set_error_temporal_data()           2020      1.5205      0.000753    1.5205      0.000753    8.34     8.34     |
|   truth_assembly()                    2000      4.2383      0.002119    4.6973      0.002349    23.26    25.77    |
|   truth_solve()                       20        1.4107      0.070535    11.3406     0.567029    7.74     62.23    |
|   update_RB_initial_condition_all_N() 20        0.0375      0.001875    0.0375      0.001875    0.21     0.21     |
|   update_RB_system_matrices()         20        0.0885      0.004425    0.3821      0.019106    0.49     2.10     |
|   update_residual_terms()             20        0.3471      0.017354    1.5892      0.079461    1.90     8.72     |
|                                                                                                                   |
| TransientRBEvaluation                                                                                             |
|   cache_online_residual_terms()       357       0.0023      0.000006    0.0023      0.000006    0.01     0.01     |
|   compute_residual_dual_norm()        35700     0.0321      0.000001    0.0321      0.000001    0.18     0.18     |
|   rb_solve()                          357       0.1013      0.000284    0.1383      0.000387    0.56     0.76     |
|   resize_data_structures()            1         0.0000      0.000045    0.0001      0.000121    0.00     0.00     |
|   write_offline_data_to_files()       1         0.0004      0.000398    0.0201      0.020124    0.00     0.11     |
 -------------------------------------------------------------------------------------------------------------------
| Totals:                               60115     18.2250                                         100.00            |
 -------------------------------------------------------------------------------------------------------------------

*** Warning, This code is untested, experimental, or likely to see future API changes: src/reduced_basis/rb_parametrized.C, line 40, compiled Aug 24 2012 at 15:15:42 ***
 EquationSystems
  n_systems()=1
   System #0, "RBConvectionDiffusion"
    Type "TransientRBConstruction"
    Variables="u" 
    Finite Element Types="LAGRANGE", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" 
    Approximation Orders="FIRST", "THIRD" 
    n_dofs()=676
    n_local_dofs()=129
    n_constrained_dofs()=102
    n_local_constrained_dofs()=21
    n_vectors()=3
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 7.94574
      Average Off-Processor Bandwidth <= 0.550388
      Maximum  On-Processor Bandwidth <= 9
      Maximum Off-Processor Bandwidth <= 5
    DofMap Constraints
      Number of DoF Constraints = 100
      Average DoF Constraint Length= 0
      Number of Node Constraints = 0

 Mesh Information:
  mesh_dimension()=2
  spatial_dimension()=3
  n_nodes()=676
    n_local_nodes()=129
  n_elem()=625
    n_local_elem()=106
    n_active_elem()=625
  n_subdomains()=1
  n_partitions()=6
  n_processors()=6
  n_threads()=1
  processor_id()=0

x_vel: 1
y_vel: 1

Error bound (absolute) at the final time is 0.0167619

************************************************************************************************************************
***             WIDEN YOUR WINDOW TO 120 CHARACTERS.  Use 'enscript -r -fCourier9' to print this document            ***
************************************************************************************************************************

---------------------------------------------- PETSc Performance Summary: ----------------------------------------------

./reduced_basis_ex3-opt on a intel-11. named daedalus with 6 processors, by roystgnr Fri Aug 24 15:22:40 2012
Using Petsc Release Version 3.1.0, Patch 5, Mon Sep 27 11:51:54 CDT 2010

                         Max       Max/Min        Avg      Total 
Time (sec):           1.135e-01      1.00942   1.127e-01
Objects:              4.600e+01      1.00000   4.600e+01
Flops:                5.160e+03      1.32990   4.507e+03  2.704e+04
Flops/sec:            4.570e+04      1.32474   3.997e+04  2.398e+05
MPI Messages:         4.000e+01      2.50000   2.533e+01  1.520e+02
MPI Message Lengths:  2.160e+03      1.96007   5.676e+01  8.628e+03
MPI Reductions:       1.260e+02      1.00000

Flop counting convention: 1 flop = 1 real number operation of type (multiply/divide/add/subtract)
                            e.g., VecAXPY() for real vectors of length N --> 2N flops
                            and VecAXPY() for complex vectors of length N --> 8N flops

Summary of Stages:   ----- Time ------  ----- Flops -----  --- Messages ---  -- Message Lengths --  -- Reductions --
                        Avg     %Total     Avg     %Total   counts   %Total     Avg         %Total   counts   %Total 
 0:      Main Stage: 1.1268e-01 100.0%  2.7040e+04 100.0%  1.520e+02 100.0%  5.676e+01      100.0%  8.900e+01  70.6% 

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

VecCopy                3 1.0 1.0991e-0457.6 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecSet                29 1.0 1.2159e-05 1.3 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecAXPY               20 1.0 4.6730e-05 1.7 5.16e+03 1.3 0.0e+00 0.0e+00 0.0e+00  0100  0  0  0   0100  0  0  0   579
VecAssemblyBegin      23 1.0 3.1102e-02 3.4 0.00e+00 0.0 0.0e+00 0.0e+00 6.9e+01 24  0  0  0 55  24  0  0  0 78     0
VecAssemblyEnd        23 1.0 5.5552e-05 2.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecScatterBegin        2 1.0 4.2200e-05 1.4 0.00e+00 0.0 3.8e+01 1.3e+02 0.0e+00  0  0 25 56  0   0  0 25 56  0     0
VecScatterEnd          2 1.0 2.4080e-0417.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatZeroEntries         2 1.0 1.0967e-05 1.4 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
------------------------------------------------------------------------------------------------------------------------

Memory usage is given in bytes:

Object Type          Creations   Destructions     Memory  Descendants' Mem.
Reports information only for process 0.

--- Event Stage 0: Main Stage

                 Vec    31             31        71016     0
         Vec Scatter     3              3         2604     0
           Index Set     6              6         3420     0
   IS L to G Mapping     3              3         3084     0
              Matrix     3              3        22784     0
========================================================================================================================
Average time to get PetscTime(): 9.53674e-08
Average time for MPI_Barrier(): 4.67777e-05
Average time for zero size MPI_Send(): 0.000114004
#PETSc Option Table entries:
-ksp_right_pc
-log_summary
-online_mode 1
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
| Time:           Fri Aug 24 15:22:40 2012                         |
| OS:             Linux                                            |
| HostName:       daedalus                                         |
| OS Release:     2.6.32-34-generic                                |
| OS Version:     #76-Ubuntu SMP Tue Aug 30 17:05:01 UTC 2011      |
| Machine:        x86_64                                           |
| Username:       roystgnr                                         |
| Configuration:  ./configure run on Wed Aug 22 12:44:06 CDT 2012  |
-------------------------------------------------------------------
 --------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=0.180623, Active time=0.077992                                               |
 --------------------------------------------------------------------------------------------------------------
| Event                            nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                            w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|--------------------------------------------------------------------------------------------------------------|
|                                                                                                              |
|                                                                                                              |
| DofMap                                                                                                       |
|   add_neighbors_to_send_list()   1         0.0001      0.000077    0.0001      0.000113    0.10     0.14     |
|   build_sparsity()               1         0.0006      0.000641    0.0009      0.000853    0.82     1.09     |
|   create_dof_constraints()       1         0.0009      0.000922    0.0012      0.001174    1.18     1.51     |
|   distribute_dofs()              1         0.0002      0.000239    0.0008      0.000840    0.31     1.08     |
|   dof_indices()                  1511      0.0005      0.000000    0.0005      0.000000    0.59     0.59     |
|   prepare_send_list()            1         0.0000      0.000003    0.0000      0.000003    0.00     0.00     |
|   reinit()                       1         0.0005      0.000454    0.0005      0.000454    0.58     0.58     |
|                                                                                                              |
| EquationSystems                                                                                              |
|   build_solution_vector()        2         0.0003      0.000157    0.0013      0.000641    0.40     1.64     |
|                                                                                                              |
| ExodusII_IO                                                                                                  |
|   write_nodal_data()             2         0.0019      0.000964    0.0019      0.000964    2.47     2.47     |
|                                                                                                              |
| Mesh                                                                                                         |
|   find_neighbors()               1         0.0005      0.000544    0.0021      0.002104    0.70     2.70     |
|   renumber_nodes_and_elem()      2         0.0000      0.000024    0.0000      0.000024    0.06     0.06     |
|                                                                                                              |
| MeshCommunication                                                                                            |
|   assign_global_indices()        1         0.0061      0.006095    0.0183      0.018281    7.81     23.44    |
|   compute_hilbert_indices()      2         0.0024      0.001223    0.0024      0.001223    3.14     3.14     |
|   find_global_indices()          2         0.0003      0.000138    0.0113      0.005634    0.35     14.45    |
|   parallel_sort()                2         0.0030      0.001483    0.0072      0.003596    3.80     9.22     |
|                                                                                                              |
| MeshOutput                                                                                                   |
|   write_equation_systems()       2         0.0000      0.000014    0.0032      0.001619    0.03     4.15     |
|                                                                                                              |
| MeshTools::Generation                                                                                        |
|   build_cube()                   1         0.0003      0.000258    0.0003      0.000258    0.33     0.33     |
|                                                                                                              |
| MetisPartitioner                                                                                             |
|   partition()                    1         0.0016      0.001571    0.0080      0.008043    2.01     10.31    |
|                                                                                                              |
| Parallel                                                                                                     |
|   allgather()                    53        0.0195      0.000368    0.0195      0.000368    25.00    25.00    |
|   barrier()                      1         0.0023      0.002259    0.0023      0.002259    2.90     2.90     |
|   broadcast()                    314       0.0006      0.000002    0.0006      0.000002    0.78     0.73     |
|   gather()                       1         0.0000      0.000045    0.0000      0.000045    0.06     0.06     |
|   max(scalar)                    2         0.0016      0.000824    0.0016      0.000824    2.11     2.11     |
|   max(vector)                    3         0.0002      0.000067    0.0002      0.000067    0.26     0.26     |
|   min(vector)                    3         0.0010      0.000349    0.0010      0.000349    1.34     1.34     |
|   probe()                        70        0.0070      0.000100    0.0070      0.000100    8.94     8.94     |
|   receive()                      70        0.0001      0.000002    0.0071      0.000102    0.18     9.12     |
|   send()                         70        0.0001      0.000001    0.0001      0.000001    0.08     0.08     |
|   send_receive()                 78        0.0002      0.000002    0.0074      0.000095    0.21     9.49     |
|   sum()                          19        0.0066      0.000347    0.0066      0.000347    8.45     8.45     |
|                                                                                                              |
| Parallel::Request                                                                                            |
|   wait()                         70        0.0000      0.000001    0.0000      0.000001    0.06     0.06     |
|                                                                                                              |
| Partitioner                                                                                                  |
|   set_node_processor_ids()       1         0.0001      0.000098    0.0057      0.005732    0.13     7.35     |
|   set_parent_processor_ids()     1         0.0000      0.000040    0.0000      0.000040    0.05     0.05     |
|                                                                                                              |
| RBConstruction                                                                                               |
|   clear()                        3         0.0002      0.000058    0.0002      0.000058    0.22     0.22     |
|   load_basis_function()          1         0.0077      0.007694    0.0077      0.007694    9.87     9.87     |
|                                                                                                              |
| RBEvaluation                                                                                                 |
|   clear()                        2         0.0000      0.000009    0.0000      0.000009    0.02     0.02     |
|   read_offline_data_from_files() 1         0.0007      0.000695    0.0008      0.000834    0.89     1.07     |
|   resize_data_structures()       1         0.0001      0.000085    0.0001      0.000085    0.11     0.11     |
|                                                                                                              |
| TransientRBConstruction                                                                                      |
|   load_rb_solution()             1         0.0003      0.000297    0.0003      0.000297    0.38     0.38     |
|                                                                                                              |
| TransientRBEvaluation                                                                                        |
|   cache_online_residual_terms()  1         0.0000      0.000022    0.0000      0.000022    0.03     0.03     |
|   compute_residual_dual_norm()   100       0.0002      0.000002    0.0002      0.000002    0.30     0.30     |
|   rb_solve()                     1         0.0095      0.009456    0.0097      0.009718    12.12    12.46    |
|   read_offline_data_from_files() 1         0.0006      0.000576    0.0014      0.001410    0.74     1.81     |
|   resize_data_structures()       1         0.0001      0.000053    0.0001      0.000138    0.07     0.18     |
 --------------------------------------------------------------------------------------------------------------
| Totals:                          2404      0.0780                                          100.00            |
 --------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running  mpirun -np 6 ./reduced_basis_ex3-opt
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
