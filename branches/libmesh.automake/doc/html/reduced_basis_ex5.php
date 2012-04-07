<?php $root=""; ?>
<?php require($root."navigation.php"); ?>
<html>
<head>
  <?php load_style($root); ?>
</head>
 
<body>
 
<?php make_navigation("reduced_basis_ex5",$root)?>
 
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
<h1>Reduced Basis: Example 5 - Reduced Cantilever Problem</h1>
Reduced Basis version of systems_of_equations_ex4: 2D cantilever


<br><br>In this example we consider the same problem as systems_of_equations_ex4,
but we introduce one parameter, which is the thickness of the cantilever.
(Note that for simplicity we do not consider a rigorous lower bound for
the coercivity constant --- to compute this bound one can use the rbOOmit
SCM classes.)

<br><br>We consider two parameters in this problem:
mu_0: scale the mesh in the y-direction
mu_1: the traction in the x-direction on the right boundary of the cantilever
mu_2: the traction in the y-direction on the right boundary of the cantilever


<br><br>Define a function to scale the mesh according to the parameter.
</div>

<div class ="fragment">
<pre>
        void scale_mesh_and_plot(EquationSystems& es, std::vector&lt;Real&gt;& mu, const std::string& filename);
        
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
        
          const unsigned int dim = 2;
        
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
          libmesh_example_assert(dim &lt;= LIBMESH_DIM, "2D support");
        
          std::string parameters_filename = "reduced_basis_ex5.in";
          GetPot infile(parameters_filename);
        
          unsigned int n_elem_x = infile("n_elem_x", 1);
          unsigned int n_elem_y = infile("n_elem_y", 1);
          Real x_size           = infile("x_size", 1.);
          Real y_size           = infile("y_size", 1.);
        
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
                                               n_elem_x, n_elem_y,
                                               0., x_size,
                                               0., y_size,
                                               QUAD9);
        
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
              equation_systems.add_system&lt;SimpleRBConstruction&gt; ("RBElasticity");
        
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
            rb_eval.rb_solve( rb_eval.get_n_basis_functions() );
        
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
              RBConstruction& rb_con = equation_systems.get_system&lt;RBConstruction&gt;("RBElasticity");
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
              scale_mesh_and_plot(equation_systems, online_mu_vector, "RB_displacement.e");
        #endif
            }
          }
        
          return 0;
        }
        
        void scale_mesh_and_plot(EquationSystems& es, std::vector&lt;Real&gt;& mu, const std::string& filename)
        {
</pre>
</div>
<div class = "comment">
Loop over the mesh nodes and move them!
</div>

<div class ="fragment">
<pre>
          MeshBase& mesh = es.get_mesh();
        
          MeshBase::const_node_iterator       node_it  = mesh.nodes_begin();
          const MeshBase::const_node_iterator node_end = mesh.nodes_end();
          
          for( ; node_it != node_end; node_it++)
          {
            Node* node = *node_it;
            
            Real y = (*node)(1);
        
            (*node)(1) = y*mu[0];
          }
        
        #ifdef LIBMESH_HAVE_EXODUS_API
          ExodusII_IO (mesh).write_equation_systems (filename, es);
        #endif
          
</pre>
</div>
<div class = "comment">
Loop over the mesh nodes and move them!
</div>

<div class ="fragment">
<pre>
          node_it  = mesh.nodes_begin();
          
          for( ; node_it != node_end; node_it++)
          {
            Node* node = *node_it;
            
            Real y = (*node)(1);
        
            (*node)(1) = 1./mu[0];
          }
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
  
  
  
  <B><FONT COLOR="#228B22">void</FONT></B> scale_mesh_and_plot(EquationSystems&amp; es, std::vector&lt;Real&gt;&amp; mu, <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp; filename);
  
  <B><FONT COLOR="#228B22">int</FONT></B> main (<B><FONT COLOR="#228B22">int</FONT></B> argc, <B><FONT COLOR="#228B22">char</FONT></B>** argv)
  {
    LibMeshInit init (argc, argv);
  
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> dim = 2;
  
  #<B><FONT COLOR="#A020F0">if</FONT></B> !defined(LIBMESH_HAVE_XDR)
    libmesh_example_assert(false, <B><FONT COLOR="#BC8F8F">&quot;--enable-xdr&quot;</FONT></B>);
  #elif defined(LIBMESH_DEFAULT_SINGLE_PRECISION)
    libmesh_example_assert(false, <B><FONT COLOR="#BC8F8F">&quot;--disable-singleprecision&quot;</FONT></B>);
  #endif
    libmesh_example_assert(libMesh::default_solver_package() == PETSC_SOLVERS, <B><FONT COLOR="#BC8F8F">&quot;--enable-petsc&quot;</FONT></B>);
  
    libmesh_example_assert(dim &lt;= LIBMESH_DIM, <B><FONT COLOR="#BC8F8F">&quot;2D support&quot;</FONT></B>);
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::string parameters_filename = <B><FONT COLOR="#BC8F8F">&quot;reduced_basis_ex5.in&quot;</FONT></B>;
    GetPot infile(parameters_filename);
  
    <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_elem_x = infile(<B><FONT COLOR="#BC8F8F">&quot;n_elem_x&quot;</FONT></B>, 1);
    <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_elem_y = infile(<B><FONT COLOR="#BC8F8F">&quot;n_elem_y&quot;</FONT></B>, 1);
    Real x_size           = infile(<B><FONT COLOR="#BC8F8F">&quot;x_size&quot;</FONT></B>, 1.);
    Real y_size           = infile(<B><FONT COLOR="#BC8F8F">&quot;y_size&quot;</FONT></B>, 1.);
  
    <B><FONT COLOR="#228B22">bool</FONT></B> store_basis_functions = infile(<B><FONT COLOR="#BC8F8F">&quot;store_basis_functions&quot;</FONT></B>, true); <I><FONT COLOR="#B22222">// Do we write the RB basis functions to disk?
</FONT></I>  
    GetPot command_line (argc, argv);
    <B><FONT COLOR="#228B22">int</FONT></B> online_mode = 0;
    <B><FONT COLOR="#A020F0">if</FONT></B> ( command_line.search(1, <B><FONT COLOR="#BC8F8F">&quot;-online_mode&quot;</FONT></B>) )
      online_mode = command_line.next(online_mode);
  
    Mesh mesh (dim);
    <B><FONT COLOR="#5F9EA0">MeshTools</FONT></B>::Generation::build_square (mesh,
                                         n_elem_x, n_elem_y,
                                         0., x_size,
                                         0., y_size,
                                         QUAD9);
  
    EquationSystems equation_systems (mesh);
  
    SimpleRBEvaluation rb_eval;
  
    <B><FONT COLOR="#A020F0">if</FONT></B>(!online_mode) <I><FONT COLOR="#B22222">// Perform the Offline stage of the RB method
</FONT></I>    {
      SimpleRBConstruction &amp; rb_con =
        equation_systems.add_system&lt;SimpleRBConstruction&gt; (<B><FONT COLOR="#BC8F8F">&quot;RBElasticity&quot;</FONT></B>);
  
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
      
      <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_parameters = infile(<B><FONT COLOR="#BC8F8F">&quot;n_parameters&quot;</FONT></B>,1);
      <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;Real&gt; online_mu_vector(n_parameters);
      <B><FONT COLOR="#A020F0">for</FONT></B>(<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i&lt;n_parameters; i++)
      {
        online_mu_vector[i] = infile(<B><FONT COLOR="#BC8F8F">&quot;online_mu&quot;</FONT></B>, online_mu_vector[i], i);
      }
  
      rb_eval.set_current_parameters(online_mu_vector);
      rb_eval.print_current_parameters();
  
      rb_eval.rb_solve( rb_eval.get_n_basis_functions() );
  
      <B><FONT COLOR="#A020F0">if</FONT></B>(store_basis_functions)
      {
        equation_systems.read(<B><FONT COLOR="#BC8F8F">&quot;equation_systems.dat&quot;</FONT></B>, READ);
        RBConstruction&amp; rb_con = equation_systems.get_system&lt;RBConstruction&gt;(<B><FONT COLOR="#BC8F8F">&quot;RBElasticity&quot;</FONT></B>);
        rb_con.rb_eval = &amp;rb_eval;
  
        rb_eval.read_in_basis_functions(rb_con);
        
        rb_con.load_rb_solution();
  #ifdef LIBMESH_HAVE_EXODUS_API
        scale_mesh_and_plot(equation_systems, online_mu_vector, <B><FONT COLOR="#BC8F8F">&quot;RB_displacement.e&quot;</FONT></B>);
  #endif
      }
    }
  
    <B><FONT COLOR="#A020F0">return</FONT></B> 0;
  }
  
  <B><FONT COLOR="#228B22">void</FONT></B> scale_mesh_and_plot(EquationSystems&amp; es, std::vector&lt;Real&gt;&amp; mu, <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp; filename)
  {
    MeshBase&amp; mesh = es.get_mesh();
  
    <B><FONT COLOR="#5F9EA0">MeshBase</FONT></B>::const_node_iterator       node_it  = mesh.nodes_begin();
    <B><FONT COLOR="#228B22">const</FONT></B> MeshBase::const_node_iterator node_end = mesh.nodes_end();
    
    <B><FONT COLOR="#A020F0">for</FONT></B>( ; node_it != node_end; node_it++)
    {
      Node* node = *node_it;
      
      Real y = (*node)(1);
  
      (*node)(1) = y*mu[0];
    }
  
  #ifdef LIBMESH_HAVE_EXODUS_API
    ExodusII_IO (mesh).write_equation_systems (filename, es);
  #endif
    
    node_it  = mesh.nodes_begin();
    
    <B><FONT COLOR="#A020F0">for</FONT></B>( ; node_it != node_end; node_it++)
    {
      Node* node = *node_it;
      
      Real y = (*node)(1);
  
      (*node)(1) = 1./mu[0];
    }
  }
</pre> 
<a name="output"></a> 
<br><br><br> <h1> The console output of the program: </h1> 
<pre>
Updated .depend
Compiling C++ (in optimized mode) reduced_basis_ex5.C...
Linking reduced_basis_ex5-opt...
***************************************************************
* Running  ./reduced_basis_ex5-opt
***************************************************************
 
*** Warning, This code is untested, experimental, or likely to see future API changes: src/reduced_basis/rb_parametrized_object.C, line 31, compiled Apr  7 2012 at 15:51:42 ***
 EquationSystems
  n_systems()=1
   System #0, "RBElasticity"
    Type "RBConstruction"
    Variables="u" "v" 
    Finite Element Types="LAGRANGE", "JACOBI_20_00" "LAGRANGE", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" "CARTESIAN" 
    Approximation Orders="SECOND", "THIRD" "SECOND", "THIRD" 
    n_dofs()=4242
    n_local_dofs()=4242
    n_constrained_dofs()=42
    n_local_constrained_dofs()=42
    n_vectors()=1
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 30.628
      Average Off-Processor Bandwidth <= 0
      Maximum  On-Processor Bandwidth <= 50
      Maximum Off-Processor Bandwidth <= 0
    DofMap Constraints
      Number of DoF Constraints = 42
      Average DoF Constraint Length= 0
      Number of Node Constraints = 0

 Mesh Information:
  mesh_dimension()=2
  spatial_dimension()=3
  n_nodes()=2121
    n_local_nodes()=2121
  n_elem()=500
    n_local_elem()=500
    n_active_elem()=500
  n_subdomains()=1
  n_partitions()=1
  n_processors()=1
  n_threads()=1
  processor_id()=0

Initializing training parameters with random training set...
Parameter 0: log scaling = 1
Parameter 1: log scaling = 0
Parameter 2: log scaling = 0


RBConstruction parameters:
system name: RBElasticity
constrained_problem: 0
Nmax: 15
Basis training error tolerance: 0.001
A_q operators attached: 3
F_q functions attached: 2
n_outputs: 0
Number of parameters: 3
Parameter 0: Min = 0.5, Max = 2
Parameter 1: Min = -1, Max = 1
Parameter 2: Min = -1, Max = 1
n_training_samples: 100
single-matrix mode? 0
reuse preconditioner? 1
use a relative error bound in greedy? 1
write out data during basis training? 0
quiet mode? 1
parameter initialized to: 
mu[0] = 0.5
mu[1] = -1
mu[2] = -1


---- Performing Greedy basis enrichment ----

---- Basis dimension: 0 ----
Performing RB solves on training set
Maximum (relative) error bound is inf

Performing truth solve at parameter:
mu[0] = 0.772358
mu[1] = 0.884244
mu[2] = 0.0699808

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 1 ----
Performing RB solves on training set
Maximum (relative) error bound is 8.95836

Performing truth solve at parameter:
mu[0] = 1.46047
mu[1] = -0.82615
mu[2] = 0.193106

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 2 ----
Performing RB solves on training set
Maximum (relative) error bound is 3.9299

Performing truth solve at parameter:
mu[0] = 1.25721
mu[1] = -0.926689
mu[2] = -0.0695853

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 3 ----
Performing RB solves on training set
Maximum (relative) error bound is 0.330827

Performing truth solve at parameter:
mu[0] = 0.690529
mu[1] = -0.781569
mu[2] = 0.039813

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 4 ----
Performing RB solves on training set
Maximum (relative) error bound is 0.0587815

Performing truth solve at parameter:
mu[0] = 0.512761
mu[1] = 0.141972
mu[2] = 0.534986

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 5 ----
Performing RB solves on training set
Maximum (relative) error bound is 0.0116996

Performing truth solve at parameter:
mu[0] = 1.99223
mu[1] = -0.399658
mu[2] = -0.242387

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 6 ----
Performing RB solves on training set
Maximum (relative) error bound is 0.00531057

Performing truth solve at parameter:
mu[0] = 1.01308
mu[1] = -0.892211
mu[2] = 0.251325

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 7 ----
Performing RB solves on training set
Maximum (relative) error bound is 0.00229261

Performing truth solve at parameter:
mu[0] = 1.96949
mu[1] = -0.975587
mu[2] = 0.865347

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 8 ----
Performing RB solves on training set
Maximum (relative) error bound is 0.00121006

Performing truth solve at parameter:
mu[0] = 0.887242
mu[1] = -0.583442
mu[2] = -0.00375969

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 9 ----
Performing RB solves on training set
Maximum (relative) error bound is 0.000214682

Specified error tolerance reached.
Writing out the basis functions...

-------------------------------------------------------------------
| Time:           Sat Apr  7 16:02:54 2012                         |
| OS:             Linux                                            |
| HostName:       lkirk-home                                       |
| OS Release:     3.0.0-17-generic                                 |
| OS Version:     #30-Ubuntu SMP Thu Mar 8 20:45:39 UTC 2012       |
| Machine:        x86_64                                           |
| Username:       benkirk                                          |
| Configuration:  ./configure run on Sat Apr  7 15:49:27 CDT 2012  |
-------------------------------------------------------------------
 --------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=18.4147, Active time=18.1771                                                 |
 --------------------------------------------------------------------------------------------------------------
| Event                            nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                            w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|--------------------------------------------------------------------------------------------------------------|
|                                                                                                              |
|                                                                                                              |
| DofMap                                                                                                       |
|   add_neighbors_to_send_list()   1         0.0005      0.000524    0.0005      0.000524    0.00     0.00     |
|   build_constraint_matrix()      3000      0.0033      0.000001    0.0033      0.000001    0.02     0.02     |
|   cnstrn_elem_mat_vec()          3000      0.0026      0.000001    0.0026      0.000001    0.01     0.01     |
|   compute_sparsity()             1         0.0156      0.015618    0.0189      0.018915    0.09     0.10     |
|   create_dof_constraints()       1         0.0044      0.004389    0.0071      0.007107    0.02     0.04     |
|   distribute_dofs()              1         0.0017      0.001724    0.0049      0.004941    0.01     0.03     |
|   dof_indices()                  10500     0.0166      0.000002    0.0166      0.000002    0.09     0.09     |
|   prepare_send_list()            1         0.0000      0.000001    0.0000      0.000001    0.00     0.00     |
|   reinit()                       1         0.0032      0.003215    0.0032      0.003215    0.02     0.02     |
|                                                                                                              |
| EquationSystems                                                                                              |
|   write()                        1         0.0217      0.021725    0.0221      0.022080    0.12     0.12     |
|                                                                                                              |
| FE                                                                                                           |
|   compute_affine_map()           3740      0.0064      0.000002    0.0064      0.000002    0.04     0.04     |
|   compute_face_map()             740       0.0045      0.000006    0.0122      0.000017    0.02     0.07     |
|   compute_shape_functions()      3740      0.0046      0.000001    0.0046      0.000001    0.03     0.03     |
|   init_face_shape_functions()    26        0.0001      0.000003    0.0001      0.000003    0.00     0.00     |
|   init_shape_functions()         746       0.0099      0.000013    0.0099      0.000013    0.05     0.05     |
|   inverse_map()                  2220      0.0074      0.000003    0.0074      0.000003    0.04     0.04     |
|                                                                                                              |
| Mesh                                                                                                         |
|   find_neighbors()               1         0.0016      0.001618    0.0016      0.001618    0.01     0.01     |
|   renumber_nodes_and_elem()      2         0.0003      0.000133    0.0003      0.000133    0.00     0.00     |
|                                                                                                              |
| MeshCommunication                                                                                            |
|   assign_global_indices()        1         0.0283      0.028350    0.0284      0.028366    0.16     0.16     |
|                                                                                                              |
| MeshTools::Generation                                                                                        |
|   build_cube()                   1         0.0022      0.002165    0.0022      0.002165    0.01     0.01     |
|                                                                                                              |
| Parallel                                                                                                     |
|   allgather()                    5         0.0000      0.000001    0.0000      0.000001    0.00     0.00     |
|   receive()                      80        0.0002      0.000003    0.0002      0.000003    0.00     0.00     |
|   send()                         80        0.0008      0.000010    0.0008      0.000010    0.00     0.00     |
|   send_receive()                 4         0.0000      0.000003    0.0000      0.000003    0.00     0.00     |
|                                                                                                              |
| Partitioner                                                                                                  |
|   single_partition()             1         0.0002      0.000191    0.0002      0.000191    0.00     0.00     |
|                                                                                                              |
| PetscLinearSolver                                                                                            |
|   solve()                        38        14.1540     0.372474    14.1540     0.372474    77.87    77.87    |
|                                                                                                              |
| RBConstruction                                                                                               |
|   add_scaled_matrix_and_vector() 6         0.1128      0.018804    0.1671      0.027843    0.62     0.92     |
|   clear()                        1         0.0003      0.000316    0.0003      0.000316    0.00     0.00     |
|   compute_Fq_representor_norms() 1         2.9638      2.963800    3.6747      3.674687    16.31    20.22    |
|   compute_max_error_bound()      10        0.0039      0.000391    0.0183      0.001835    0.02     0.10     |
|   enrich_RB_space()              9         0.0336      0.003735    0.0336      0.003735    0.18     0.18     |
|   train_reduced_basis()          1         0.0027      0.002730    17.9265     17.926470   0.02     98.62    |
|   truth_assembly()               9         0.1226      0.013625    0.1226      0.013625    0.67     0.67     |
|   truth_solve()                  9         0.0092      0.001028    3.9224      0.435824    0.05     21.58    |
|   update_RB_system_matrices()    9         0.2090      0.023226    0.2090      0.023226    1.15     1.15     |
|   update_residual_terms()        9         0.4130      0.045885    10.0656     1.118395    2.27     55.37    |
|                                                                                                              |
| RBEvaluation                                                                                                 |
|   clear()                        1         0.0000      0.000049    0.0001      0.000112    0.00     0.00     |
|   clear_riesz_representors()     2         0.0001      0.000032    0.0001      0.000032    0.00     0.00     |
|   compute_residual_dual_norm()   1000      0.0081      0.000008    0.0081      0.000008    0.04     0.04     |
|   rb_solve()                     1000      0.0060      0.000006    0.0143      0.000014    0.03     0.08     |
|   resize_data_structures()       1         0.0001      0.000058    0.0001      0.000059    0.00     0.00     |
|   write_offline_data_to_files()  1         0.0017      0.001674    0.0017      0.001674    0.01     0.01     |
 --------------------------------------------------------------------------------------------------------------
| Totals:                          30001     18.1771                                         100.00            |
 --------------------------------------------------------------------------------------------------------------

*** Warning, This code is untested, experimental, or likely to see future API changes: src/reduced_basis/rb_parametrized_object.C, line 31, compiled Apr  7 2012 at 15:51:42 ***
mu[0] = 0.5
mu[1] = 0.5
mu[2] = -0.5

Reading in the basis functions...
Finished reading in the basis functions...

-------------------------------------------------------------------
| Time:           Sat Apr  7 16:02:54 2012                         |
| OS:             Linux                                            |
| HostName:       lkirk-home                                       |
| OS Release:     3.0.0-17-generic                                 |
| OS Version:     #30-Ubuntu SMP Thu Mar 8 20:45:39 UTC 2012       |
| Machine:        x86_64                                           |
| Username:       benkirk                                          |
| Configuration:  ./configure run on Sat Apr  7 15:49:27 CDT 2012  |
-------------------------------------------------------------------
 --------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=0.281891, Active time=0.07565                                                |
 --------------------------------------------------------------------------------------------------------------
| Event                            nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                            w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|--------------------------------------------------------------------------------------------------------------|
|                                                                                                              |
|                                                                                                              |
| DofMap                                                                                                       |
|   add_neighbors_to_send_list()   1         0.0005      0.000497    0.0005      0.000497    0.66     0.66     |
|   compute_sparsity()             1         0.0068      0.006785    0.0076      0.007560    8.97     9.99     |
|   create_dof_constraints()       1         0.0003      0.000301    0.0003      0.000301    0.40     0.40     |
|   distribute_dofs()              1         0.0016      0.001638    0.0048      0.004811    2.17     6.36     |
|   dof_indices()                  1500      0.0013      0.000001    0.0013      0.000001    1.72     1.72     |
|   prepare_send_list()            1         0.0000      0.000001    0.0000      0.000001    0.00     0.00     |
|   reinit()                       1         0.0032      0.003170    0.0032      0.003170    4.19     4.19     |
|                                                                                                              |
| EquationSystems                                                                                              |
|   build_solution_vector()        1         0.0013      0.001345    0.0021      0.002131    1.78     2.82     |
|   read()                         1         0.0245      0.024525    0.0654      0.065368    32.42    86.41    |
|   update()                       1         0.0001      0.000066    0.0001      0.000066    0.09     0.09     |
|                                                                                                              |
| ExodusII_IO                                                                                                  |
|   write_nodal_data()             1         0.0021      0.002091    0.0021      0.002091    2.76     2.76     |
|                                                                                                              |
| Mesh                                                                                                         |
|   find_neighbors()               1         0.0015      0.001532    0.0015      0.001532    2.03     2.03     |
|   renumber_nodes_and_elem()      2         0.0002      0.000122    0.0002      0.000122    0.32     0.32     |
|                                                                                                              |
| MeshCommunication                                                                                            |
|   assign_global_indices()        1         0.0276      0.027595    0.0276      0.027610    36.48    36.50    |
|                                                                                                              |
| MeshOutput                                                                                                   |
|   write_equation_systems()       1         0.0000      0.000037    0.0043      0.004260    0.05     5.63     |
|                                                                                                              |
| MeshTools::Generation                                                                                        |
|   build_cube()                   1         0.0014      0.001389    0.0014      0.001389    1.84     1.84     |
|                                                                                                              |
| Parallel                                                                                                     |
|   allgather()                    45        0.0000      0.000001    0.0000      0.000001    0.06     0.06     |
|   send_receive()                 4         0.0000      0.000002    0.0000      0.000002    0.01     0.01     |
|                                                                                                              |
| Partitioner                                                                                                  |
|   single_partition()             1         0.0002      0.000169    0.0002      0.000169    0.22     0.22     |
|                                                                                                              |
| RBConstruction                                                                                               |
|   clear()                        2         0.0002      0.000088    0.0002      0.000088    0.23     0.23     |
|   load_rb_solution()             1         0.0002      0.000180    0.0002      0.000180    0.24     0.24     |
|                                                                                                              |
| RBEvaluation                                                                                                 |
|   clear()                        1         0.0001      0.000061    0.0001      0.000062    0.08     0.08     |
|   clear_riesz_representors()     2         0.0000      0.000003    0.0000      0.000003    0.01     0.01     |
|   compute_residual_dual_norm()   1         0.0000      0.000046    0.0000      0.000046    0.06     0.06     |
|   rb_solve()                     1         0.0005      0.000504    0.0006      0.000551    0.67     0.73     |
|   read_offline_data_from_files() 1         0.0018      0.001808    0.0019      0.001938    2.39     2.56     |
|   resize_data_structures()       1         0.0001      0.000124    0.0001      0.000130    0.16     0.17     |
 --------------------------------------------------------------------------------------------------------------
| Totals:                          1576      0.0756                                          100.00            |
 --------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running  ./reduced_basis_ex5-opt
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
