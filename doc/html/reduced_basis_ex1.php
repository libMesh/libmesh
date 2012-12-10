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
We need to give the RBConstruction object a pointer to
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
    
  <I><FONT COLOR="#B22222">/* rbOOmit is distributed in the hope that it will be useful, */</FONT></I>
  <I><FONT COLOR="#B22222">/* but WITHOUT ANY WARRANTY; without even the implied warranty of */</FONT></I>
  <I><FONT COLOR="#B22222">/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU */</FONT></I>
  <I><FONT COLOR="#B22222">/* Lesser General Public License for more details. */</FONT></I>
    
  <I><FONT COLOR="#B22222">/* You should have received a copy of the GNU Lesser General Public */</FONT></I>
  <I><FONT COLOR="#B22222">/* License along with this library; if not, write to the Free Software */</FONT></I>
  <I><FONT COLOR="#B22222">/* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA */</FONT></I>
  
   
  
  
  
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
Linking reduced_basis_ex1-opt...
***************************************************************
* Running  mpirun -np 6 ./reduced_basis_ex1-opt -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
***************************************************************
 
*** Warning, This code is untested, experimental, or likely to see future API changes: src/reduced_basis/rb_parametrized.C, line 40, compiled Aug 24 2012 at 15:15:42 ***
 EquationSystems
  n_systems()=1
   System #0, "RBConvectionDiffusion"
    Type "RBConstruction"
    Variables="u" 
    Finite Element Types="LAGRANGE", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" 
    Approximation Orders="FIRST", "THIRD" 
    n_dofs()=2601
    n_local_dofs()=462
    n_constrained_dofs()=203
    n_local_constrained_dofs()=42
    n_vectors()=1
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 8.45022
      Average Off-Processor Bandwidth <= 0.279221
      Maximum  On-Processor Bandwidth <= 9
      Maximum Off-Processor Bandwidth <= 5
    DofMap Constraints
      Number of DoF Constraints = 200
      Average DoF Constraint Length= 0
      Number of Node Constraints = 0

 Mesh Information:
  mesh_dimension()=2
  spatial_dimension()=3
  n_nodes()=2601
    n_local_nodes()=462
  n_elem()=2500
    n_local_elem()=417
    n_active_elem()=2500
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
Parameter x_vel: Min = -2, Max = 2, value = -2
Parameter y_vel: Min = -2, Max = 2, value = 0.79
n_training_samples: 100
single-matrix mode? 0
reuse preconditioner? 1
use a relative error bound in greedy? 0
write out data during basis training? 0
quiet mode? 1

Compute output dual norms
output_dual_innerprods[0][0] = 0.323675
output_dual_innerprods[1][0] = 0.323675
output_dual_innerprods[2][0] = 0.323675
output_dual_innerprods[3][0] = 0.323675

---- Performing Greedy basis enrichment ----

---- Basis dimension: 0 ----
Performing RB solves on training set
Maximum (absolute) error bound is 3.74824

Performing truth solve at parameter:
x_vel: -2
y_vel: -2

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 1 ----
Performing RB solves on training set
Maximum (absolute) error bound is 6.67625

Performing truth solve at parameter:
x_vel: 2
y_vel: 2

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 2 ----
Performing RB solves on training set
Maximum (absolute) error bound is 5.18295

Performing truth solve at parameter:
x_vel: -2
y_vel: 2

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 3 ----
Performing RB solves on training set
Maximum (absolute) error bound is 3.59028

Performing truth solve at parameter:
x_vel: 2
y_vel: -2

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 4 ----
Performing RB solves on training set
Maximum (absolute) error bound is 2.71338

Performing truth solve at parameter:
x_vel: -0.222222
y_vel: 0.222222

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 5 ----
Performing RB solves on training set
Maximum (absolute) error bound is 2.50784

Performing truth solve at parameter:
x_vel: 0.666667
y_vel: -0.222222

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 6 ----
Performing RB solves on training set
Maximum (absolute) error bound is 2.15653

Performing truth solve at parameter:
x_vel: -0.222222
y_vel: -0.666667

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 7 ----
Performing RB solves on training set
Maximum (absolute) error bound is 1.71352

Performing truth solve at parameter:
x_vel: 0.222222
y_vel: 1.11111

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 8 ----
Performing RB solves on training set
Maximum (absolute) error bound is 1.61442

Performing truth solve at parameter:
x_vel: -1.55556
y_vel: -0.222222

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 9 ----
Performing RB solves on training set
Maximum (absolute) error bound is 1.09712

Performing truth solve at parameter:
x_vel: 2
y_vel: 0.222222

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 10 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.9794

Performing truth solve at parameter:
x_vel: 0.222222
y_vel: -2

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 11 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.763958

Performing truth solve at parameter:
x_vel: 0.222222
y_vel: 0.222222

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 12 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.616746

Performing truth solve at parameter:
x_vel: -1.11111
y_vel: 1.11111

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 13 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.530181

Performing truth solve at parameter:
x_vel: -0.222222
y_vel: 2

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 14 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.475932

Performing truth solve at parameter:
x_vel: -0.666667
y_vel: -0.222222

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 15 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.440399

Performing truth solve at parameter:
x_vel: 1.11111
y_vel: 0.666667

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 16 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.423537

Performing truth solve at parameter:
x_vel: -1.11111
y_vel: -1.11111

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 17 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.408558

Performing truth solve at parameter:
x_vel: 0.666667
y_vel: -1.11111

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 18 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.333764

Performing truth solve at parameter:
x_vel: 0.222222
y_vel: -0.222222

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 19 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.196867

Performing truth solve at parameter:
x_vel: -0.222222
y_vel: 0.666667

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 20 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.191475

Maximum number of basis functions reached: Nmax = 20
************************************************************************************************************************
***             WIDEN YOUR WINDOW TO 120 CHARACTERS.  Use 'enscript -r -fCourier9' to print this document            ***
************************************************************************************************************************

---------------------------------------------- PETSc Performance Summary: ----------------------------------------------

./reduced_basis_ex1-opt on a intel-11. named daedalus with 6 processors, by roystgnr Fri Aug 24 15:22:08 2012
Using Petsc Release Version 3.1.0, Patch 5, Mon Sep 27 11:51:54 CDT 2010

                         Max       Max/Min        Avg      Total 
Time (sec):           2.872e+00      1.00025   2.872e+00
Objects:              4.960e+02      1.20388   4.400e+02
Flops:                4.145e+08      1.18394   3.853e+08  2.312e+09
Flops/sec:            1.443e+08      1.18394   1.342e+08  8.051e+08
MPI Messages:         5.422e+04      2.50000   3.254e+04  1.952e+05
MPI Message Lengths:  9.816e+06      2.36863   1.770e+02  3.456e+07
MPI Reductions:       3.252e+04      1.00259

Flop counting convention: 1 flop = 1 real number operation of type (multiply/divide/add/subtract)
                            e.g., VecAXPY() for real vectors of length N --> 2N flops
                            and VecAXPY() for complex vectors of length N --> 8N flops

Summary of Stages:   ----- Time ------  ----- Flops -----  --- Messages ---  -- Message Lengths --  -- Reductions --
                        Avg     %Total     Avg     %Total   counts   %Total     Avg         %Total   counts   %Total 
 0:      Main Stage: 2.8715e+00 100.0%  2.3118e+09 100.0%  1.952e+05 100.0%  1.770e+02      100.0%  3.218e+04  99.2% 

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

VecDot              4075 1.0 3.0059e-01 1.9 3.76e+06 1.2 0.0e+00 0.0e+00 4.1e+03  8  1  0  0 13   8  1  0  0 13    70
VecMDot             6580 1.0 6.8266e-01 1.8 8.83e+07 1.2 0.0e+00 0.0e+00 6.6e+03 17 22  0  0 20  17 22  0  0 20   729
VecNorm             6907 1.0 4.7556e-01 1.3 6.38e+06 1.2 0.0e+00 0.0e+00 6.9e+03 15  2  0  0 21  15  2  0  0 21    76
VecScale            6922 1.0 4.2653e-03 1.2 3.19e+06 1.2 0.0e+00 0.0e+00 0.0e+00  0  1  0  0  0   0  1  0  0  0  4209
VecCopy             1106 1.0 5.7983e-04 1.3 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecSet             12507 1.0 5.8296e-03 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecAXPY              699 1.0 6.3562e-04 1.2 6.46e+05 1.2 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0  5721
VecMAXPY            6822 1.0 3.7299e-02 1.2 9.45e+07 1.2 0.0e+00 0.0e+00 0.0e+00  1 23  0  0  0   1 23  0  0  0 14267
VecAssemblyBegin     595 1.0 1.1685e-01 1.1 0.00e+00 0.0 9.0e+01 2.8e+02 1.8e+03  4  0  0  0  6   4  0  0  0  6     0
VecAssemblyEnd       595 1.0 3.5572e-04 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecScatterBegin    10822 1.0 3.4482e-02 1.7 0.00e+00 0.0 1.9e+05 1.8e+02 0.0e+00  1  0100100  0   1  0100100  0     0
VecScatterEnd      10822 1.0 4.9496e-0113.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00 10  0  0  0  0  10  0  0  0  0     0
VecNormalize        6822 1.0 4.7867e-01 1.3 9.46e+06 1.2 0.0e+00 0.0e+00 6.8e+03 15  2  0  0 21  15  2  0  0 21   111
MatMult             6822 1.0 4.4858e-01 4.2 5.19e+07 1.2 1.2e+05 1.8e+02 0.0e+00 10 13 63 63  0  10 13 63 63  0   654
MatMultAdd          3915 1.0 1.7491e-01 4.5 3.16e+07 1.2 7.0e+04 1.8e+02 0.0e+00  4  8 36 36  0   4  8 36 36  0  1021
MatSolve            6822 1.0 9.5048e-02 1.2 1.30e+08 1.2 0.0e+00 0.0e+00 0.0e+00  3 31  0  0  0   3 31  0  0  0  7480
MatLUFactorNum        42 1.0 7.5040e-03 1.2 3.85e+06 1.3 0.0e+00 0.0e+00 0.0e+00  0  1  0  0  0   0  1  0  0  0  2767
MatILUFactorSym       42 1.0 1.7970e-02 1.2 0.00e+00 0.0 0.0e+00 0.0e+00 4.2e+01  1  0  0  0  0   1  0  0  0  0     0
MatAssemblyBegin    4214 1.0 3.9824e-01 1.1 0.00e+00 0.0 1.1e+02 9.0e+02 8.4e+03 13  0  0  0 26  13  0  0  0 26     0
MatAssemblyEnd      4214 1.0 2.5899e-01 1.2 0.00e+00 0.0 1.8e+02 4.6e+01 4.2e+03  8  0  0  0 13   8  0  0  0 13     0
MatGetRow          37884 1.2 4.6742e-03 1.3 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatGetRowIJ           42 1.0 2.2173e-05 3.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatGetOrdering        42 1.0 9.7680e-04 2.7 0.00e+00 0.0 0.0e+00 0.0e+00 1.1e+02  0  0  0  0  0   0  0  0  0  0     0
MatZeroEntries        56 1.0 2.6870e-04 1.6 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPGMRESOrthog      6580 1.0 7.2229e-01 1.7 1.77e+08 1.2 0.0e+00 0.0e+00 6.6e+03 18 43  0  0 20  18 43  0  0 20  1378
KSPSetup             127 1.0 9.7752e-05 1.5 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPSolve              85 1.0 1.4203e+00 1.0 3.79e+08 1.2 1.2e+05 1.8e+02 1.4e+04 49 91 63 63 42  49 91 63 63 42  1486
PCSetUp               84 1.0 2.7427e-02 1.2 3.85e+06 1.3 0.0e+00 0.0e+00 1.5e+02  1  1  0  0  0   1  1  0  0  0   757
PCSetUpOnBlocks       85 1.0 2.7054e-02 1.2 3.85e+06 1.3 0.0e+00 0.0e+00 1.5e+02  1  1  0  0  0   1  1  0  0  0   767
PCApply             6822 1.0 1.4142e-01 1.2 1.30e+08 1.2 0.0e+00 0.0e+00 0.0e+00  5 31  0  0  0   5 31  0  0  0  5027
------------------------------------------------------------------------------------------------------------------------

Memory usage is given in bytes:

Object Type          Creations   Destructions     Memory  Descendants' Mem.
Reports information only for process 0.

--- Event Stage 0: Main Stage

                 Vec   202            202       957392     0
         Vec Scatter     7              7         6076     0
           Index Set   140            140       306864     0
   IS L to G Mapping     1              1         2456     0
              Matrix    58             58      5384476     0
       Krylov Solver     2              2        18880     0
      Preconditioner     2              2         1408     0
========================================================================================================================
Average time to get PetscTime(): 9.53674e-08
Average time for MPI_Barrier(): 4.4632e-05
Average time for zero size MPI_Send(): 4.56572e-05
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
| Time:           Fri Aug 24 15:22:08 2012                         |
| OS:             Linux                                            |
| HostName:       daedalus                                         |
| OS Release:     2.6.32-34-generic                                |
| OS Version:     #76-Ubuntu SMP Tue Aug 30 17:05:01 UTC 2011      |
| Machine:        x86_64                                           |
| Username:       roystgnr                                         |
| Configuration:  ./configure run on Wed Aug 22 12:44:06 CDT 2012  |
-------------------------------------------------------------------
 -------------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=3.06119, Active time=2.84514                                                      |
 -------------------------------------------------------------------------------------------------------------------
| Event                                 nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                                 w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|-------------------------------------------------------------------------------------------------------------------|
|                                                                                                                   |
|                                                                                                                   |
| DofMap                                                                                                            |
|   add_neighbors_to_send_list()        1         0.0002      0.000230    0.0003      0.000260    0.01     0.01     |
|   build_constraint_matrix()           3753      0.0018      0.000000    0.0018      0.000000    0.06     0.06     |
|   build_sparsity()                    1         0.0024      0.002420    0.0033      0.003254    0.09     0.11     |
|   cnstrn_elem_mat_vec()               3753      0.0082      0.000002    0.0082      0.000002    0.29     0.29     |
|   create_dof_constraints()            1         0.0029      0.002924    0.0037      0.003728    0.10     0.13     |
|   distribute_dofs()                   1         0.0006      0.000643    0.0051      0.005131    0.02     0.18     |
|   dof_indices()                       12603     0.0043      0.000000    0.0043      0.000000    0.15     0.15     |
|   prepare_send_list()                 1         0.0000      0.000006    0.0000      0.000006    0.00     0.00     |
|   reinit()                            1         0.0015      0.001463    0.0015      0.001463    0.05     0.05     |
|                                                                                                                   |
| FE                                                                                                                |
|   compute_shape_functions()           8226      0.0187      0.000002    0.0187      0.000002    0.66     0.66     |
|   init_shape_functions()              738       0.0019      0.000003    0.0019      0.000003    0.07     0.07     |
|   inverse_map()                       1440      0.0014      0.000001    0.0014      0.000001    0.05     0.05     |
|                                                                                                                   |
| FEMap                                                                                                             |
|   compute_affine_map()                8226      0.0049      0.000001    0.0049      0.000001    0.17     0.17     |
|   compute_face_map()                  720       0.0017      0.000002    0.0032      0.000004    0.06     0.11     |
|   init_face_shape_functions()         18        0.0000      0.000002    0.0000      0.000002    0.00     0.00     |
|   init_reference_to_physical_map()    738       0.0011      0.000001    0.0011      0.000001    0.04     0.04     |
|                                                                                                                   |
| Mesh                                                                                                              |
|   find_neighbors()                    1         0.0019      0.001945    0.0056      0.005612    0.07     0.20     |
|   renumber_nodes_and_elem()           2         0.0002      0.000079    0.0002      0.000079    0.01     0.01     |
|                                                                                                                   |
| MeshCommunication                                                                                                 |
|   assign_global_indices()             1         0.0203      0.020336    0.0364      0.036442    0.71     1.28     |
|   compute_hilbert_indices()           2         0.0088      0.004398    0.0088      0.004398    0.31     0.31     |
|   find_global_indices()               2         0.0008      0.000382    0.0215      0.010774    0.03     0.76     |
|   parallel_sort()                     2         0.0013      0.000668    0.0103      0.005166    0.05     0.36     |
|                                                                                                                   |
| MeshTools::Generation                                                                                             |
|   build_cube()                        1         0.0007      0.000712    0.0007      0.000712    0.03     0.03     |
|                                                                                                                   |
| MetisPartitioner                                                                                                  |
|   partition()                         1         0.0040      0.003971    0.0145      0.014523    0.14     0.51     |
|                                                                                                                   |
| Parallel                                                                                                          |
|   allgather()                         13        0.0167      0.001287    0.0167      0.001287    0.59     0.59     |
|   barrier()                           21        0.0014      0.000065    0.0014      0.000065    0.05     0.05     |
|   broadcast()                         43        0.0008      0.000019    0.0008      0.000019    0.03     0.03     |
|   gather()                            81        0.0083      0.000103    0.0083      0.000103    0.29     0.29     |
|   max(scalar)                         2         0.0081      0.004026    0.0081      0.004026    0.28     0.28     |
|   max(vector)                         3         0.0001      0.000042    0.0001      0.000042    0.00     0.00     |
|   maxloc(scalar)                      21        0.0413      0.001965    0.0413      0.001965    1.45     1.45     |
|   min(vector)                         3         0.0008      0.000283    0.0008      0.000283    0.03     0.03     |
|   probe()                             70        0.0065      0.000093    0.0065      0.000093    0.23     0.23     |
|   receive()                           550       0.0004      0.000001    0.0069      0.000013    0.01     0.24     |
|   send()                              150       0.0002      0.000001    0.0002      0.000001    0.01     0.01     |
|   send_receive()                      78        0.0001      0.000002    0.0069      0.000089    0.01     0.24     |
|   sum()                               34        0.0123      0.000363    0.0123      0.000363    0.43     0.43     |
|                                                                                                                   |
| Parallel::Request                                                                                                 |
|   wait()                              550       0.0005      0.000001    0.0005      0.000001    0.02     0.02     |
|                                                                                                                   |
| Partitioner                                                                                                       |
|   set_node_processor_ids()            1         0.0003      0.000275    0.0048      0.004814    0.01     0.17     |
|   set_parent_processor_ids()          1         0.0001      0.000146    0.0001      0.000146    0.01     0.01     |
|                                                                                                                   |
| PetscLinearSolver                                                                                                 |
|   solve()                             85        1.4684      0.017276    1.4684      0.017276    51.61    51.61    |
|                                                                                                                   |
| RBConstruction                                                                                                    |
|   add_scaled_matrix_and_vector()      9         0.0824      0.009157    0.1278      0.014196    2.90     4.49     |
|   clear()                             1         0.0005      0.000482    0.0005      0.000482    0.02     0.02     |
|   compute_Fq_representor_innerprods() 1         0.0100      0.010005    0.0200      0.019968    0.35     0.70     |
|   compute_max_error_bound()           21        0.0013      0.000060    0.1007      0.004794    0.04     3.54     |
|   compute_output_dual_innerprods()    1         0.0131      0.013114    0.0879      0.087933    0.46     3.09     |
|   enrich_RB_space()                   20        0.0795      0.003975    0.0795      0.003975    2.79     2.79     |
|   train_reduced_basis()               1         0.0154      0.015415    2.5817      2.581723    0.54     90.74    |
|   truth_assembly()                    20        0.0402      0.002009    0.0402      0.002009    1.41     1.41     |
|   truth_solve()                       20        0.0143      0.000714    0.3361      0.016806    0.50     11.81    |
|   update_RB_system_matrices()         20        0.2850      0.014250    0.2850      0.014250    10.02    10.02    |
|   update_residual_terms()             20        0.5550      0.027751    1.6570      0.082851    19.51    58.24    |
|                                                                                                                   |
| RBEvaluation                                                                                                      |
|   clear()                             1         0.0001      0.000083    0.0001      0.000083    0.00     0.00     |
|   compute_residual_dual_norm()        357       0.0501      0.000140    0.0501      0.000140    1.76     1.76     |
|   rb_solve()                          357       0.0057      0.000016    0.0558      0.000156    0.20     1.96     |
|   resize_data_structures()            1         0.0001      0.000068    0.0001      0.000068    0.00     0.00     |
|   write_offline_data_to_files()       1         0.0364      0.036394    0.0364      0.036394    1.28     1.28     |
 -------------------------------------------------------------------------------------------------------------------
| Totals:                               42790     2.8451                                          100.00            |
 -------------------------------------------------------------------------------------------------------------------

*** Warning, This code is untested, experimental, or likely to see future API changes: src/reduced_basis/rb_parametrized.C, line 40, compiled Aug 24 2012 at 15:15:42 ***
 EquationSystems
  n_systems()=1
   System #0, "RBConvectionDiffusion"
    Type "RBConstruction"
    Variables="u" 
    Finite Element Types="LAGRANGE", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" 
    Approximation Orders="FIRST", "THIRD" 
    n_dofs()=2601
    n_local_dofs()=462
    n_constrained_dofs()=203
    n_local_constrained_dofs()=42
    n_vectors()=1
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 8.45022
      Average Off-Processor Bandwidth <= 0.279221
      Maximum  On-Processor Bandwidth <= 9
      Maximum Off-Processor Bandwidth <= 5
    DofMap Constraints
      Number of DoF Constraints = 200
      Average DoF Constraint Length= 0
      Number of Node Constraints = 0

 Mesh Information:
  mesh_dimension()=2
  spatial_dimension()=3
  n_nodes()=2601
    n_local_nodes()=462
  n_elem()=2500
    n_local_elem()=417
    n_active_elem()=2500
  n_subdomains()=1
  n_partitions()=6
  n_processors()=6
  n_threads()=1
  processor_id()=0

x_vel: -2
y_vel: 0.79

output 1, value = 0.115628, bound = 0.104928
output 2, value = 0.371774, bound = 0.104928
output 3, value = 0.253449, bound = 0.104928
output 4, value = 0.12462, bound = 0.104928

************************************************************************************************************************
***             WIDEN YOUR WINDOW TO 120 CHARACTERS.  Use 'enscript -r -fCourier9' to print this document            ***
************************************************************************************************************************

---------------------------------------------- PETSc Performance Summary: ----------------------------------------------

./reduced_basis_ex1-opt on a intel-11. named daedalus with 6 processors, by roystgnr Fri Aug 24 15:22:09 2012
Using Petsc Release Version 3.1.0, Patch 5, Mon Sep 27 11:51:54 CDT 2010

                         Max       Max/Min        Avg      Total 
Time (sec):           1.573e-01      1.00754   1.563e-01
Objects:              3.200e+01      1.00000   3.200e+01
Flops:                1.848e+04      1.17857   1.734e+04  1.040e+05
Flops/sec:            1.175e+05      1.16977   1.109e+05  6.655e+05
MPI Messages:         2.000e+01      2.50000   1.200e+01  7.200e+01
MPI Message Lengths:  3.400e+03      2.36439   1.677e+02  1.207e+04
MPI Reductions:       1.080e+02      1.00000

Flop counting convention: 1 flop = 1 real number operation of type (multiply/divide/add/subtract)
                            e.g., VecAXPY() for real vectors of length N --> 2N flops
                            and VecAXPY() for complex vectors of length N --> 8N flops

Summary of Stages:   ----- Time ------  ----- Flops -----  --- Messages ---  -- Message Lengths --  -- Reductions --
                        Avg     %Total     Avg     %Total   counts   %Total     Avg         %Total   counts   %Total 
 0:      Main Stage: 1.5627e-01 100.0%  1.0404e+05 100.0%  7.200e+01 100.0%  1.677e+02      100.0%  7.700e+01  71.3% 

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

VecCopy                3 1.0 1.0014e-05 5.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecSet                25 1.0 1.8358e-05 2.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecAXPY               20 1.0 6.4135e-05 1.9 1.85e+04 1.2 0.0e+00 0.0e+00 0.0e+00  0100  0  0  0   0100  0  0  0  1622
VecAssemblyBegin      23 1.0 7.3469e-03 1.7 0.00e+00 0.0 0.0e+00 0.0e+00 6.9e+01  4  0  0  0 64   4  0  0  0 90     0
VecAssemblyEnd        23 1.0 4.3154e-05 1.7 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecScatterBegin        2 1.0 4.7922e-05 1.8 0.00e+00 0.0 3.6e+01 2.7e+02 0.0e+00  0  0 50 80  0   0  0 50 80  0     0
VecScatterEnd          2 1.0 3.5787e-0410.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatZeroEntries         2 1.0 1.2875e-05 1.4 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
------------------------------------------------------------------------------------------------------------------------

Memory usage is given in bytes:

Object Type          Creations   Destructions     Memory  Descendants' Mem.
Reports information only for process 0.

--- Event Stage 0: Main Stage

                 Vec    25             25       122648     0
         Vec Scatter     1              1          868     0
           Index Set     2              2         1236     0
   IS L to G Mapping     1              1         2456     0
              Matrix     3              3        66020     0
========================================================================================================================
Average time to get PetscTime(): 9.53674e-08
Average time for MPI_Barrier(): 1.55926e-05
Average time for zero size MPI_Send(): 1.01725e-05
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
| Time:           Fri Aug 24 15:22:09 2012                         |
| OS:             Linux                                            |
| HostName:       daedalus                                         |
| OS Release:     2.6.32-34-generic                                |
| OS Version:     #76-Ubuntu SMP Tue Aug 30 17:05:01 UTC 2011      |
| Machine:        x86_64                                           |
| Username:       roystgnr                                         |
| Configuration:  ./configure run on Wed Aug 22 12:44:06 CDT 2012  |
-------------------------------------------------------------------
 --------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=0.296775, Active time=0.13068                                                |
 --------------------------------------------------------------------------------------------------------------
| Event                            nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                            w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|--------------------------------------------------------------------------------------------------------------|
|                                                                                                              |
|                                                                                                              |
| DofMap                                                                                                       |
|   add_neighbors_to_send_list()   1         0.0002      0.000227    0.0003      0.000263    0.17     0.20     |
|   build_sparsity()               1         0.0024      0.002429    0.0033      0.003282    1.86     2.51     |
|   create_dof_constraints()       1         0.0028      0.002809    0.0036      0.003649    2.15     2.79     |
|   distribute_dofs()              1         0.0006      0.000641    0.0058      0.005764    0.49     4.41     |
|   dof_indices()                  5931      0.0016      0.000000    0.0016      0.000000    1.23     1.23     |
|   prepare_send_list()            1         0.0000      0.000005    0.0000      0.000005    0.00     0.00     |
|   reinit()                       1         0.0015      0.001460    0.0015      0.001460    1.12     1.12     |
|                                                                                                              |
| EquationSystems                                                                                              |
|   build_solution_vector()        2         0.0009      0.000464    0.0030      0.001519    0.71     2.32     |
|                                                                                                              |
| ExodusII_IO                                                                                                  |
|   write_nodal_data()             2         0.0038      0.001892    0.0038      0.001892    2.90     2.90     |
|                                                                                                              |
| Mesh                                                                                                         |
|   find_neighbors()               1         0.0019      0.001925    0.0058      0.005765    1.47     4.41     |
|   renumber_nodes_and_elem()      2         0.0002      0.000086    0.0002      0.000086    0.13     0.13     |
|                                                                                                              |
| MeshCommunication                                                                                            |
|   assign_global_indices()        1         0.0201      0.020099    0.0365      0.036520    15.38    27.95    |
|   compute_hilbert_indices()      2         0.0088      0.004391    0.0088      0.004391    6.72     6.72     |
|   find_global_indices()          2         0.0008      0.000387    0.0232      0.011617    0.59     17.78    |
|   parallel_sort()                2         0.0015      0.000773    0.0104      0.005200    1.18     7.96     |
|                                                                                                              |
| MeshOutput                                                                                                   |
|   write_equation_systems()       2         0.0000      0.000012    0.0068      0.003423    0.02     5.24     |
|                                                                                                              |
| MeshTools::Generation                                                                                        |
|   build_cube()                   1         0.0007      0.000731    0.0007      0.000731    0.56     0.56     |
|                                                                                                              |
| MetisPartitioner                                                                                             |
|   partition()                    1         0.0040      0.003953    0.0156      0.015580    3.02     11.92    |
|                                                                                                              |
| Parallel                                                                                                     |
|   allgather()                    53        0.0285      0.000538    0.0285      0.000538    21.82    21.82    |
|   barrier()                      1         0.0058      0.005771    0.0058      0.005771    4.42     4.42     |
|   broadcast()                    314       0.0024      0.000008    0.0024      0.000008    1.83     1.81     |
|   gather()                       1         0.0000      0.000006    0.0000      0.000006    0.00     0.00     |
|   max(scalar)                    2         0.0074      0.003707    0.0074      0.003707    5.67     5.67     |
|   max(vector)                    3         0.0001      0.000035    0.0001      0.000035    0.08     0.08     |
|   min(vector)                    3         0.0026      0.000860    0.0026      0.000860    1.97     1.97     |
|   probe()                        70        0.0078      0.000112    0.0078      0.000112    6.00     6.00     |
|   receive()                      70        0.0002      0.000002    0.0080      0.000114    0.12     6.12     |
|   send()                         70        0.0001      0.000001    0.0001      0.000001    0.06     0.06     |
|   send_receive()                 78        0.0001      0.000002    0.0083      0.000106    0.11     6.34     |
|   sum()                          19        0.0125      0.000658    0.0125      0.000658    9.57     9.57     |
|                                                                                                              |
| Parallel::Request                                                                                            |
|   wait()                         70        0.0000      0.000001    0.0000      0.000001    0.04     0.04     |
|                                                                                                              |
| Partitioner                                                                                                  |
|   set_node_processor_ids()       1         0.0003      0.000269    0.0035      0.003455    0.21     2.64     |
|   set_parent_processor_ids()     1         0.0001      0.000145    0.0001      0.000145    0.11     0.11     |
|                                                                                                              |
| RBConstruction                                                                                               |
|   clear()                        1         0.0002      0.000175    0.0002      0.000175    0.13     0.13     |
|   load_basis_function()          1         0.0002      0.000189    0.0002      0.000189    0.14     0.14     |
|   load_rb_solution()             1         0.0004      0.000422    0.0004      0.000422    0.32     0.32     |
|                                                                                                              |
| RBEvaluation                                                                                                 |
|   clear()                        1         0.0001      0.000064    0.0001      0.000064    0.05     0.05     |
|   compute_residual_dual_norm()   1         0.0005      0.000460    0.0005      0.000460    0.35     0.35     |
|   rb_solve()                     1         0.0088      0.008778    0.0092      0.009238    6.72     7.07     |
|   read_offline_data_from_files() 1         0.0007      0.000671    0.0007      0.000743    0.51     0.57     |
|   resize_data_structures()       1         0.0001      0.000072    0.0001      0.000072    0.06     0.06     |
 --------------------------------------------------------------------------------------------------------------
| Totals:                          6719      0.1307                                          100.00            |
 --------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running  mpirun -np 6 ./reduced_basis_ex1-opt -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
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
