<?php $root=""; ?>
<?php require($root."navigation.php"); ?>
<html>
<head>
  <?php load_style($root); ?>
</head>
 
<body>
 
<?php make_navigation("reduced_basis_ex4",$root)?>
 
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
<h1>Reduced Basis Example 4 - Empirical Interpolation Method</h1>


<br><br>In this example problem we develop a reduced basis approximation for a parametrized
PDE that has "non-affine" parameter dependence. This requires the use of the
Empirical Interpolation Method (EIM).

<br><br>We first use EIM to construct an affine approximation to the non-affine term,
which is a parametrized function that is a Gaussian with "center" defined
by the two parameters (mu_1,mu_2) \in [-1,1]^2. We then employ this EIM
approximation in order to generate a reduced basis approximation for the
parametrized PDE: -0.05 * Laplacian(u) = f(mu_1,mu_2), with zero Dirichlet
boundary conditions.


<br><br>Basic include file needed for the mesh functionality.
</div>

<div class ="fragment">
<pre>
        #include "libmesh.h"
        #include "mesh.h"
        #include "mesh_generation.h"
        #include "equation_systems.h"
        #include "exodusII_io.h"
        #include "getpot.h"
        
        #include "eim_classes.h"
        #include "rb_classes.h"
        
        
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
Skip this 2D example if libMesh was compiled as 1D-only.
</div>

<div class ="fragment">
<pre>
          libmesh_example_assert(2 &lt;= LIBMESH_DIM, "2D support");
        
</pre>
</div>
<div class = "comment">
Define the names of the input files we will read the problem properties from
</div>

<div class ="fragment">
<pre>
          std::string eim_parameters = "eim.in";
          std::string rb_parameters  = "rb.in";
          std::string main_parameters = "reduced_basis_ex4.in";
          GetPot infile(main_parameters);
        
          unsigned int n_elem = infile("n_elem", 1);       // Determines the number of elements in the "truth" mesh
          const unsigned int dim = 2;                      // The number of spatial dimensions
          bool store_basis_functions = infile("store_basis_functions", false); // Do we write out basis functions?
        
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
Create a mesh (just a simple square)
</div>

<div class ="fragment">
<pre>
          Mesh mesh (dim);
          MeshTools::Generation::build_square (mesh,
                                               n_elem, n_elem,
                                               -1., 1.,
                                               -1., 1.,
                                               QUAD4);
        
</pre>
</div>
<div class = "comment">
Initialize the EquationSystems object for this mesh and attach
the EIM and RB Construction objects
</div>

<div class ="fragment">
<pre>
          EquationSystems equation_systems (mesh);
        
          SimpleEIMConstruction & eim_construction =
            equation_systems.add_system&lt;SimpleEIMConstruction&gt; ("EIM");
          SimpleRBConstruction & rb_construction =
            equation_systems.add_system&lt;SimpleRBConstruction&gt; ("RB");
          
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
          mesh.print_info();
          equation_systems.print_info();
        
</pre>
</div>
<div class = "comment">
Initialize the standard RBEvaluation object
</div>

<div class ="fragment">
<pre>
          SimpleRBEvaluation rb_eval;
        
</pre>
</div>
<div class = "comment">
Initialize the EIM RBEvaluation object
</div>

<div class ="fragment">
<pre>
          SimpleEIMEvaluation eim_rb_eval;
          
</pre>
</div>
<div class = "comment">
Set the rb_eval objects for the RBConstructions
</div>

<div class ="fragment">
<pre>
          eim_construction.set_rb_evaluation(eim_rb_eval);
          rb_construction.set_rb_evaluation(rb_eval);
        
          if(!online_mode)
          {
</pre>
</div>
<div class = "comment">
Read data from input file and print state
</div>

<div class ="fragment">
<pre>
            eim_construction.process_parameters_file(eim_parameters);
            eim_construction.print_info();
          
</pre>
</div>
<div class = "comment">
Perform the EIM Greedy and write out the data
</div>

<div class ="fragment">
<pre>
            eim_construction.initialize_rb_construction();
            eim_construction.train_reduced_basis();
            eim_construction.get_rb_evaluation().write_offline_data_to_files("eim_data");
        
</pre>
</div>
<div class = "comment">
Read data from input file and print state
</div>

<div class ="fragment">
<pre>
            rb_construction.process_parameters_file(rb_parameters);
        
</pre>
</div>
<div class = "comment">
attach the EIM theta objects to the RBConstruction and RBEvaluation objects
</div>

<div class ="fragment">
<pre>
            eim_rb_eval.initialize_eim_theta_objects();
            rb_eval.get_rb_theta_expansion().attach_multiple_F_theta(eim_rb_eval.get_eim_theta_objects());
            
</pre>
</div>
<div class = "comment">
attach the EIM assembly objects to the RBConstruction object
</div>

<div class ="fragment">
<pre>
            eim_construction.initialize_eim_assembly_objects();
            rb_construction.get_rb_assembly_expansion().attach_multiple_F_assembly(eim_construction.get_eim_assembly_objects());
        
</pre>
</div>
<div class = "comment">
Print out the state of rb_construction now that the EIM objects have been attached
</div>

<div class ="fragment">
<pre>
            rb_construction.print_info();
        
</pre>
</div>
<div class = "comment">
Need to initialize _after_ EIM greedy so that
the system knows how many affine terms there are
</div>

<div class ="fragment">
<pre>
            rb_construction.initialize_rb_construction();
            rb_construction.train_reduced_basis();
            rb_construction.get_rb_evaluation().write_offline_data_to_files("rb_data");
        
</pre>
</div>
<div class = "comment">
Write out the basis functions, if requested
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
              eim_construction.get_rb_evaluation().write_out_basis_functions(eim_construction,"eim_data");
              rb_construction.get_rb_evaluation().write_out_basis_functions(rb_construction,"rb_data");
            }
          }
          else
          {
            eim_rb_eval.read_offline_data_from_files("eim_data");
        
</pre>
</div>
<div class = "comment">
attach the EIM theta objects to rb_eval objects
</div>

<div class ="fragment">
<pre>
            eim_rb_eval.initialize_eim_theta_objects();
            rb_eval.get_rb_theta_expansion().attach_multiple_F_theta(eim_rb_eval.get_eim_theta_objects());
            
</pre>
</div>
<div class = "comment">
Read in the offline data for rb_eval
</div>

<div class ="fragment">
<pre>
            rb_eval.read_offline_data_from_files("rb_data");
        
</pre>
</div>
<div class = "comment">
Get the parameters at which we will do a reduced basis solve
</div>

<div class ="fragment">
<pre>
            Real online_center_x = infile("online_center_x", 0.);
            Real online_center_y = infile("online_center_y", 0.);
            RBParameters online_mu;
            online_mu.set_value("center_x", online_center_x);
            online_mu.set_value("center_y", online_center_y);
            rb_eval.set_parameters(online_mu);
            rb_eval.print_parameters();
            rb_eval.rb_solve( rb_eval.get_n_basis_functions() );
        
</pre>
</div>
<div class = "comment">
plot the solution, if requested
</div>

<div class ="fragment">
<pre>
            if(store_basis_functions)
            {
</pre>
</div>
<div class = "comment">
read in the data from files
</div>

<div class ="fragment">
<pre>
              eim_rb_eval.read_in_basis_functions(eim_construction,"eim_data");
              rb_eval.read_in_basis_functions(rb_construction,"rb_data");
        
              eim_construction.load_rb_solution();
              rb_construction.load_rb_solution();
        #ifdef LIBMESH_HAVE_EXODUS_API
              ExodusII_IO(mesh).write_equation_systems("RB_sol.e",equation_systems);
        #endif
            }
          }
        
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
  #include <B><FONT COLOR="#BC8F8F">&quot;equation_systems.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;exodusII_io.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;getpot.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;eim_classes.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;rb_classes.h&quot;</FONT></B>
  
  
  <B><FONT COLOR="#228B22">int</FONT></B> main (<B><FONT COLOR="#228B22">int</FONT></B> argc, <B><FONT COLOR="#228B22">char</FONT></B>** argv)
  {
    LibMeshInit init (argc, argv);
  
  #<B><FONT COLOR="#A020F0">if</FONT></B> !defined(LIBMESH_HAVE_XDR)
    libmesh_example_assert(false, <B><FONT COLOR="#BC8F8F">&quot;--enable-xdr&quot;</FONT></B>);
  #elif defined(LIBMESH_DEFAULT_SINGLE_PRECISION)
    libmesh_example_assert(false, <B><FONT COLOR="#BC8F8F">&quot;--disable-singleprecision&quot;</FONT></B>);
  #endif
  
    libmesh_example_assert(2 &lt;= LIBMESH_DIM, <B><FONT COLOR="#BC8F8F">&quot;2D support&quot;</FONT></B>);
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::string eim_parameters = <B><FONT COLOR="#BC8F8F">&quot;eim.in&quot;</FONT></B>;
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::string rb_parameters  = <B><FONT COLOR="#BC8F8F">&quot;rb.in&quot;</FONT></B>;
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::string main_parameters = <B><FONT COLOR="#BC8F8F">&quot;reduced_basis_ex4.in&quot;</FONT></B>;
    GetPot infile(main_parameters);
  
    <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_elem = infile(<B><FONT COLOR="#BC8F8F">&quot;n_elem&quot;</FONT></B>, 1);       <I><FONT COLOR="#B22222">// Determines the number of elements in the &quot;truth&quot; mesh
</FONT></I>    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> dim = 2;                      <I><FONT COLOR="#B22222">// The number of spatial dimensions
</FONT></I>    <B><FONT COLOR="#228B22">bool</FONT></B> store_basis_functions = infile(<B><FONT COLOR="#BC8F8F">&quot;store_basis_functions&quot;</FONT></B>, false); <I><FONT COLOR="#B22222">// Do we write out basis functions?
</FONT></I>  
    GetPot command_line (argc, argv);
    <B><FONT COLOR="#228B22">int</FONT></B> online_mode = 0;
    <B><FONT COLOR="#A020F0">if</FONT></B> ( command_line.search(1, <B><FONT COLOR="#BC8F8F">&quot;-online_mode&quot;</FONT></B>) )
      online_mode = command_line.next(online_mode);
  
    Mesh mesh (dim);
    <B><FONT COLOR="#5F9EA0">MeshTools</FONT></B>::Generation::build_square (mesh,
                                         n_elem, n_elem,
                                         -1., 1.,
                                         -1., 1.,
                                         QUAD4);
  
    EquationSystems equation_systems (mesh);
  
    SimpleEIMConstruction &amp; eim_construction =
      equation_systems.add_system&lt;SimpleEIMConstruction&gt; (<B><FONT COLOR="#BC8F8F">&quot;EIM&quot;</FONT></B>);
    SimpleRBConstruction &amp; rb_construction =
      equation_systems.add_system&lt;SimpleRBConstruction&gt; (<B><FONT COLOR="#BC8F8F">&quot;RB&quot;</FONT></B>);
    
    equation_systems.init ();
  
    mesh.print_info();
    equation_systems.print_info();
  
    SimpleRBEvaluation rb_eval;
  
    SimpleEIMEvaluation eim_rb_eval;
    
    eim_construction.set_rb_evaluation(eim_rb_eval);
    rb_construction.set_rb_evaluation(rb_eval);
  
    <B><FONT COLOR="#A020F0">if</FONT></B>(!online_mode)
    {
      eim_construction.process_parameters_file(eim_parameters);
      eim_construction.print_info();
    
      eim_construction.initialize_rb_construction();
      eim_construction.train_reduced_basis();
      eim_construction.get_rb_evaluation().write_offline_data_to_files(<B><FONT COLOR="#BC8F8F">&quot;eim_data&quot;</FONT></B>);
  
      rb_construction.process_parameters_file(rb_parameters);
  
      eim_rb_eval.initialize_eim_theta_objects();
      rb_eval.get_rb_theta_expansion().attach_multiple_F_theta(eim_rb_eval.get_eim_theta_objects());
      
      eim_construction.initialize_eim_assembly_objects();
      rb_construction.get_rb_assembly_expansion().attach_multiple_F_assembly(eim_construction.get_eim_assembly_objects());
  
      rb_construction.print_info();
  
      rb_construction.initialize_rb_construction();
      rb_construction.train_reduced_basis();
      rb_construction.get_rb_evaluation().write_offline_data_to_files(<B><FONT COLOR="#BC8F8F">&quot;rb_data&quot;</FONT></B>);
  
      <B><FONT COLOR="#A020F0">if</FONT></B>(store_basis_functions)
      {
        eim_construction.get_rb_evaluation().write_out_basis_functions(eim_construction,<B><FONT COLOR="#BC8F8F">&quot;eim_data&quot;</FONT></B>);
        rb_construction.get_rb_evaluation().write_out_basis_functions(rb_construction,<B><FONT COLOR="#BC8F8F">&quot;rb_data&quot;</FONT></B>);
      }
    }
    <B><FONT COLOR="#A020F0">else</FONT></B>
    {
      eim_rb_eval.read_offline_data_from_files(<B><FONT COLOR="#BC8F8F">&quot;eim_data&quot;</FONT></B>);
  
      eim_rb_eval.initialize_eim_theta_objects();
      rb_eval.get_rb_theta_expansion().attach_multiple_F_theta(eim_rb_eval.get_eim_theta_objects());
      
      rb_eval.read_offline_data_from_files(<B><FONT COLOR="#BC8F8F">&quot;rb_data&quot;</FONT></B>);
  
      Real online_center_x = infile(<B><FONT COLOR="#BC8F8F">&quot;online_center_x&quot;</FONT></B>, 0.);
      Real online_center_y = infile(<B><FONT COLOR="#BC8F8F">&quot;online_center_y&quot;</FONT></B>, 0.);
      RBParameters online_mu;
      online_mu.set_value(<B><FONT COLOR="#BC8F8F">&quot;center_x&quot;</FONT></B>, online_center_x);
      online_mu.set_value(<B><FONT COLOR="#BC8F8F">&quot;center_y&quot;</FONT></B>, online_center_y);
      rb_eval.set_parameters(online_mu);
      rb_eval.print_parameters();
      rb_eval.rb_solve( rb_eval.get_n_basis_functions() );
  
      <B><FONT COLOR="#A020F0">if</FONT></B>(store_basis_functions)
      {
        eim_rb_eval.read_in_basis_functions(eim_construction,<B><FONT COLOR="#BC8F8F">&quot;eim_data&quot;</FONT></B>);
        rb_eval.read_in_basis_functions(rb_construction,<B><FONT COLOR="#BC8F8F">&quot;rb_data&quot;</FONT></B>);
  
        eim_construction.load_rb_solution();
        rb_construction.load_rb_solution();
  #ifdef LIBMESH_HAVE_EXODUS_API
        ExodusII_IO(mesh).write_equation_systems(<B><FONT COLOR="#BC8F8F">&quot;RB_sol.e&quot;</FONT></B>,equation_systems);
  #endif
      }
    }
  
  }
</pre> 
<a name="output"></a> 
<br><br><br> <h1> The console output of the program: </h1> 
<pre>
Linking reduced_basis_ex4-opt...
***************************************************************
* Running  ./reduced_basis_ex4-opt
***************************************************************
 
*** Warning, This code is untested, experimental, or likely to see future API changes: src/reduced_basis/rb_parametrized.C, line 40, compiled Aug 24 2012 at 15:15:42 ***
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

 EquationSystems
  n_systems()=2
   System #0, "EIM"
    Type "RBConstruction"
    Variables="f_EIM" 
    Finite Element Types="LAGRANGE", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" 
    Approximation Orders="FIRST", "THIRD" 
    n_dofs()=676
    n_local_dofs()=676
    n_constrained_dofs()=0
    n_local_constrained_dofs()=0
    n_vectors()=1
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 8.54438
      Average Off-Processor Bandwidth <= 0
      Maximum  On-Processor Bandwidth <= 9
      Maximum Off-Processor Bandwidth <= 0
    DofMap Constraints
      Number of DoF Constraints = 0
      Number of Node Constraints = 0
   System #1, "RB"
    Type "RBConstruction"
    Variables="u" 
    Finite Element Types="LAGRANGE", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" 
    Approximation Orders="FIRST", "THIRD" 
    n_dofs()=676
    n_local_dofs()=676
    n_constrained_dofs()=100
    n_local_constrained_dofs()=100
    n_vectors()=1
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

Initializing training parameters with deterministic training set...
Parameter center_x: log scaling = 0
Parameter center_y: log scaling = 0


RBConstruction parameters:
system name: EIM
constrained_problem: 0
Nmax: 20
Basis training error tolerance: 0.001
Aq operators attached: 0
Fq functions attached: 0
n_outputs: 0
Number of parameters: 2
Parameter center_x: Min = -1, Max = 1, value = 1
Parameter center_y: Min = -1, Max = 1, value = 1
n_training_samples: 25
single-matrix mode? 0
reuse preconditioner? 1
use a relative error bound in greedy? 0
write out data during basis training? 0
quiet mode? 1


RBEIMConstruction parameters:
best fit type: projection

Initializing parametrized functions in training set...
Parametrized functions in training set initialized


---- Performing Greedy basis enrichment ----

---- Basis dimension: 0 ----
Performing truth solve at parameter:
center_x: 1
center_y: 1

Enriching the RB space
Updating RB matrices

---- Basis dimension: 1 ----
Performing RB solves on training set
Maximum (absolute) error bound is 1.04028

Performing truth solve at parameter:
center_x: 0.5
center_y: 0.5

Enriching the RB space
Updating RB matrices

---- Basis dimension: 2 ----
Performing RB solves on training set
Maximum (absolute) error bound is 1.00428

Performing truth solve at parameter:
center_x: -1
center_y: -1

Enriching the RB space
Updating RB matrices

---- Basis dimension: 3 ----
Performing RB solves on training set
Maximum (absolute) error bound is 1.02817

Performing truth solve at parameter:
center_x: -0.5
center_y: -0.5

Enriching the RB space
Updating RB matrices

---- Basis dimension: 4 ----
Performing RB solves on training set
Maximum (absolute) error bound is 1.00298

Performing truth solve at parameter:
center_x: 1
center_y: -1

Enriching the RB space
Updating RB matrices

---- Basis dimension: 5 ----
Performing RB solves on training set
Maximum (absolute) error bound is 1.00292

Performing truth solve at parameter:
center_x: -1
center_y: 1

Enriching the RB space
Updating RB matrices

---- Basis dimension: 6 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.86594

Performing truth solve at parameter:
center_x: 0.5
center_y: -0.5

Enriching the RB space
Updating RB matrices

---- Basis dimension: 7 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.816903

Performing truth solve at parameter:
center_x: -0.5
center_y: 0.5

Enriching the RB space
Updating RB matrices

---- Basis dimension: 8 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.681622

Performing truth solve at parameter:
center_x: 0
center_y: -1

Enriching the RB space
Updating RB matrices

---- Basis dimension: 9 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.658046

Performing truth solve at parameter:
center_x: 0
center_y: 1

Enriching the RB space
Updating RB matrices

---- Basis dimension: 10 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.606155

Performing truth solve at parameter:
center_x: 1
center_y: 0

Enriching the RB space
Updating RB matrices

---- Basis dimension: 11 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.606138

Performing truth solve at parameter:
center_x: -1
center_y: 0

Enriching the RB space
Updating RB matrices

---- Basis dimension: 12 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.193607

Performing truth solve at parameter:
center_x: 0
center_y: 0

Enriching the RB space
Updating RB matrices

---- Basis dimension: 13 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.184868

Performing truth solve at parameter:
center_x: -0.5
center_y: 1

Enriching the RB space
Updating RB matrices

---- Basis dimension: 14 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.214039

Performing truth solve at parameter:
center_x: -1
center_y: 0.5

Enriching the RB space
Updating RB matrices

---- Basis dimension: 15 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.175349

Performing truth solve at parameter:
center_x: 1
center_y: 0.5

Enriching the RB space
Updating RB matrices

---- Basis dimension: 16 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.163258

Performing truth solve at parameter:
center_x: 0.5
center_y: 1

Enriching the RB space
Updating RB matrices

---- Basis dimension: 17 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.161335

Performing truth solve at parameter:
center_x: 0.5
center_y: -1

Enriching the RB space
Updating RB matrices

---- Basis dimension: 18 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.159752

Performing truth solve at parameter:
center_x: 1
center_y: -0.5

Enriching the RB space
Updating RB matrices

---- Basis dimension: 19 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.123948

Performing truth solve at parameter:
center_x: -0.5
center_y: -1

Enriching the RB space
Updating RB matrices

---- Basis dimension: 20 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.162888

Maximum number of basis functions reached: Nmax = 20.
Perform one more Greedy iteration for error bounds.
Performing truth solve at parameter:
center_x: -1
center_y: -0.5

Enriching the RB space
Updating RB matrices

---- Basis dimension: 20 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.162888

Extra Greedy iteration finished.
Initializing training parameters with deterministic training set...
Parameter center_x: log scaling = 0
Parameter center_y: log scaling = 0


RBConstruction parameters:
system name: RB
constrained_problem: 0
Nmax: 15
Basis training error tolerance: 0.001
Aq operators attached: 1
Fq functions attached: 20
n_outputs: 0
Number of parameters: 2
Parameter center_x: Min = -1, Max = 1, value = 1
Parameter center_y: Min = -1, Max = 1, value = 1
n_training_samples: 100
single-matrix mode? 0
reuse preconditioner? 1
use a relative error bound in greedy? 0
write out data during basis training? 0
quiet mode? 1


---- Performing Greedy basis enrichment ----

---- Basis dimension: 0 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.388596

Performing truth solve at parameter:
center_x: 0.111111
center_y: -0.111111

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 1 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.176922

Performing truth solve at parameter:
center_x: -0.555556
center_y: 0.555556

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 2 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.156429

Performing truth solve at parameter:
center_x: 0.555556
center_y: 0.555556

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 3 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.130432

Performing truth solve at parameter:
center_x: -0.555556
center_y: -0.555556

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 4 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.0995698

Performing truth solve at parameter:
center_x: 0.555556
center_y: -0.555556

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 5 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.0476481

Performing truth solve at parameter:
center_x: 1
center_y: 0.111111

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 6 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.0433854

Performing truth solve at parameter:
center_x: -1
center_y: -0.111111

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 7 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.0309096

Performing truth solve at parameter:
center_x: 0.111111
center_y: -1

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 8 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.0266356

Performing truth solve at parameter:
center_x: -0.111111
center_y: 1

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 9 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.0243385

Performing truth solve at parameter:
center_x: -1
center_y: 1

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 10 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.0230817

Performing truth solve at parameter:
center_x: -1
center_y: -1

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 11 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.020751

Performing truth solve at parameter:
center_x: 1
center_y: 1

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 12 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.0181426

Performing truth solve at parameter:
center_x: 1
center_y: -1

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 13 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.00823617

Performing truth solve at parameter:
center_x: -0.333333
center_y: -1

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 14 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.00764269

Performing truth solve at parameter:
center_x: 0.333333
center_y: 1

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 15 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.00557789

Maximum number of basis functions reached: Nmax = 15

-------------------------------------------------------------------
| Time:           Fri Aug 24 15:22:45 2012                         |
| OS:             Linux                                            |
| HostName:       daedalus                                         |
| OS Release:     2.6.32-34-generic                                |
| OS Version:     #76-Ubuntu SMP Tue Aug 30 17:05:01 UTC 2011      |
| Machine:        x86_64                                           |
| Username:       roystgnr                                         |
| Configuration:  ./configure run on Wed Aug 22 12:44:06 CDT 2012  |
-------------------------------------------------------------------
 -------------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=1.8024, Active time=1.64069                                                       |
 -------------------------------------------------------------------------------------------------------------------
| Event                                 nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                                 w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|-------------------------------------------------------------------------------------------------------------------|
|                                                                                                                   |
|                                                                                                                   |
| DofMap                                                                                                            |
|   add_neighbors_to_send_list()        2         0.0002      0.000088    0.0002      0.000088    0.01     0.01     |
|   build_constraint_matrix()           13750     0.0078      0.000001    0.0078      0.000001    0.47     0.47     |
|   build_sparsity()                    2         0.0016      0.000802    0.0021      0.001030    0.10     0.13     |
|   cnstrn_elem_mat_vec()               13750     0.0076      0.000001    0.0076      0.000001    0.46     0.46     |
|   create_dof_constraints()            2         0.0009      0.000450    0.0011      0.000556    0.05     0.07     |
|   distribute_dofs()                   2         0.0003      0.000142    0.0011      0.000536    0.02     0.07     |
|   dof_indices()                       101065    0.0340      0.000000    0.0340      0.000000    2.07     2.07     |
|   prepare_send_list()                 2         0.0000      0.000000    0.0000      0.000000    0.00     0.00     |
|   reinit()                            2         0.0008      0.000392    0.0008      0.000392    0.05     0.05     |
|                                                                                                                   |
| FE                                                                                                                |
|   compute_shape_functions()           90850     0.1143      0.000001    0.1143      0.000001    6.96     6.96     |
|   init_shape_functions()              4738      0.0108      0.000002    0.0108      0.000002    0.66     0.66     |
|   inverse_map()                       60080     0.0502      0.000001    0.0502      0.000001    3.06     3.06     |
|                                                                                                                   |
| FEMap                                                                                                             |
|   compute_affine_map()                90850     0.0504      0.000001    0.0504      0.000001    3.07     3.07     |
|   compute_face_map()                  4600      0.0099      0.000002    0.0184      0.000004    0.60     1.12     |
|   init_face_shape_functions()         46        0.0001      0.000002    0.0001      0.000002    0.01     0.01     |
|   init_reference_to_physical_map()    4738      0.0067      0.000001    0.0067      0.000001    0.41     0.41     |
|                                                                                                                   |
| Mesh                                                                                                              |
|   find_neighbors()                    1         0.0005      0.000505    0.0005      0.000505    0.03     0.03     |
|   renumber_nodes_and_elem()           2         0.0000      0.000019    0.0000      0.000019    0.00     0.00     |
|                                                                                                                   |
| MeshCommunication                                                                                                 |
|   assign_global_indices()             2         0.0139      0.006974    0.0140      0.006983    0.85     0.85     |
|                                                                                                                   |
| MeshTools::Generation                                                                                             |
|   build_cube()                        1         0.0002      0.000238    0.0002      0.000238    0.01     0.01     |
|                                                                                                                   |
| Parallel                                                                                                          |
|   allgather()                         10        0.0000      0.000001    0.0000      0.000001    0.00     0.00     |
|   receive()                           140       0.0002      0.000001    0.0002      0.000001    0.01     0.01     |
|   send()                              140       0.0002      0.000002    0.0002      0.000002    0.01     0.01     |
|   send_receive()                      8         0.0000      0.000001    0.0000      0.000001    0.00     0.00     |
|                                                                                                                   |
| Parallel::Request                                                                                                 |
|   wait()                              140       0.0001      0.000001    0.0001      0.000001    0.01     0.01     |
|                                                                                                                   |
| Partitioner                                                                                                       |
|   single_partition()                  1         0.0000      0.000026    0.0000      0.000026    0.00     0.00     |
|                                                                                                                   |
| PetscLinearSolver                                                                                                 |
|   solve()                             75        0.0825      0.001100    0.0825      0.001100    5.03     5.03     |
|                                                                                                                   |
| PointLocatorTree                                                                                                  |
|   init(no master)                     1         0.0004      0.000400    0.0004      0.000400    0.02     0.02     |
|   operator()                          440       0.0034      0.000008    0.0039      0.000009    0.20     0.24     |
|                                                                                                                   |
| RBConstruction                                                                                                    |
|   add_scaled_matrix_and_vector()      23        0.1501      0.006525    0.3346      0.014550    9.15     20.40    |
|   clear()                             3         0.0004      0.000126    0.0004      0.000126    0.02     0.02     |
|   compute_Fq_representor_innerprods() 2         0.0056      0.002812    0.0311      0.015546    0.34     1.90     |
|   compute_max_error_bound()           37        0.0066      0.000179    0.8308      0.022453    0.40     50.64    |
|   enrich_RB_space()                   15        0.0027      0.000177    0.0027      0.000177    0.16     0.16     |
|   train_reduced_basis()               2         0.0019      0.000928    1.0492      0.524609    0.11     63.95    |
|   truth_assembly()                    15        0.0037      0.000249    0.0041      0.000273    0.23     0.25     |
|   truth_solve()                       15        0.0005      0.000033    0.0260      0.001735    0.03     1.59     |
|   update_RB_system_matrices()         36        0.0084      0.000232    0.0084      0.000232    0.51     0.51     |
|   update_residual_terms()             15        0.0105      0.000697    0.0324      0.002158    0.64     1.97     |
|                                                                                                                   |
| RBEIMConstruction                                                                                                 |
|   compute_best_fit_error()            525       0.1070      0.000204    0.1139      0.000217    6.52     6.94     |
|   enrich_RB_space()                   21        0.0495      0.002357    0.1120      0.005331    3.02     6.82     |
|   truth_solve()                       571       0.1308      0.000229    0.2177      0.000381    7.97     13.27    |
|   update_RB_system_matrices()         21        0.0024      0.000115    0.0080      0.000380    0.15     0.49     |
|                                                                                                                   |
| RBEIMEvaluation                                                                                                   |
|   rb_solve()                          1633      0.0286      0.000017    0.0286      0.000017    1.74     1.74     |
|   write_offline_data_to_files()       1         0.0002      0.000203    0.0273      0.027278    0.01     1.66     |
|                                                                                                                   |
| RBEvaluation                                                                                                      |
|   clear()                             3         0.0000      0.000011    0.0000      0.000011    0.00     0.00     |
|   compute_residual_dual_norm()        1600      0.6460      0.000404    0.6460      0.000404    39.37    39.37    |
|   rb_solve()                          1600      0.0358      0.000022    0.7100      0.000444    2.18     43.28    |
|   resize_data_structures()            2         0.0001      0.000028    0.0001      0.000028    0.00     0.00     |
|   write_offline_data_to_files()       2         0.0530      0.026509    0.0530      0.026509    3.23     3.23     |
 -------------------------------------------------------------------------------------------------------------------
| Totals:                               391582    1.6407                                          100.00            |
 -------------------------------------------------------------------------------------------------------------------

*** Warning, This code is untested, experimental, or likely to see future API changes: src/reduced_basis/rb_parametrized.C, line 40, compiled Aug 24 2012 at 15:15:42 ***
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

 EquationSystems
  n_systems()=2
   System #0, "EIM"
    Type "RBConstruction"
    Variables="f_EIM" 
    Finite Element Types="LAGRANGE", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" 
    Approximation Orders="FIRST", "THIRD" 
    n_dofs()=676
    n_local_dofs()=676
    n_constrained_dofs()=0
    n_local_constrained_dofs()=0
    n_vectors()=1
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 8.54438
      Average Off-Processor Bandwidth <= 0
      Maximum  On-Processor Bandwidth <= 9
      Maximum Off-Processor Bandwidth <= 0
    DofMap Constraints
      Number of DoF Constraints = 0
      Number of Node Constraints = 0
   System #1, "RB"
    Type "RBConstruction"
    Variables="u" 
    Finite Element Types="LAGRANGE", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" 
    Approximation Orders="FIRST", "THIRD" 
    n_dofs()=676
    n_local_dofs()=676
    n_constrained_dofs()=100
    n_local_constrained_dofs()=100
    n_vectors()=1
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

center_x: -0.6
center_y: 0.7


-------------------------------------------------------------------
| Time:           Fri Aug 24 15:22:45 2012                         |
| OS:             Linux                                            |
| HostName:       daedalus                                         |
| OS Release:     2.6.32-34-generic                                |
| OS Version:     #76-Ubuntu SMP Tue Aug 30 17:05:01 UTC 2011      |
| Machine:        x86_64                                           |
| Username:       roystgnr                                         |
| Configuration:  ./configure run on Wed Aug 22 12:44:06 CDT 2012  |
-------------------------------------------------------------------
 --------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=0.048012, Active time=0.030355                                               |
 --------------------------------------------------------------------------------------------------------------
| Event                            nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                            w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|--------------------------------------------------------------------------------------------------------------|
|                                                                                                              |
|                                                                                                              |
| DofMap                                                                                                       |
|   add_neighbors_to_send_list()   2         0.0002      0.000088    0.0002      0.000088    0.58     0.58     |
|   build_sparsity()               2         0.0016      0.000817    0.0021      0.001034    5.39     6.81     |
|   create_dof_constraints()       2         0.0009      0.000455    0.0011      0.000561    3.00     3.70     |
|   distribute_dofs()              2         0.0003      0.000138    0.0011      0.000532    0.91     3.51     |
|   dof_indices()                  3125      0.0008      0.000000    0.0008      0.000000    2.79     2.79     |
|   prepare_send_list()            2         0.0000      0.000000    0.0000      0.000000    0.00     0.00     |
|   reinit()                       2         0.0008      0.000392    0.0008      0.000392    2.59     2.59     |
|                                                                                                              |
| EquationSystems                                                                                              |
|   build_solution_vector()        1         0.0010      0.000995    0.0014      0.001415    3.28     4.66     |
|                                                                                                              |
| ExodusII_IO                                                                                                  |
|   write_nodal_data()             1         0.0013      0.001290    0.0013      0.001290    4.25     4.25     |
|                                                                                                              |
| Mesh                                                                                                         |
|   find_neighbors()               1         0.0005      0.000502    0.0005      0.000502    1.65     1.65     |
|   renumber_nodes_and_elem()      2         0.0000      0.000019    0.0000      0.000019    0.13     0.13     |
|                                                                                                              |
| MeshCommunication                                                                                            |
|   assign_global_indices()        2         0.0138      0.006923    0.0139      0.006934    45.62    45.68    |
|                                                                                                              |
| MeshOutput                                                                                                   |
|   write_equation_systems()       1         0.0000      0.000016    0.0027      0.002721    0.05     8.96     |
|                                                                                                              |
| MeshTools::Generation                                                                                        |
|   build_cube()                   1         0.0002      0.000235    0.0002      0.000235    0.77     0.77     |
|                                                                                                              |
| Parallel                                                                                                     |
|   allgather()                    80        0.0000      0.000000    0.0000      0.000000    0.09     0.09     |
|   send_receive()                 8         0.0000      0.000002    0.0000      0.000002    0.05     0.05     |
|                                                                                                              |
| Partitioner                                                                                                  |
|   single_partition()             1         0.0000      0.000025    0.0000      0.000025    0.08     0.08     |
|                                                                                                              |
| RBConstruction                                                                                               |
|   clear()                        3         0.0002      0.000063    0.0002      0.000063    0.62     0.62     |
|   load_rb_solution()             2         0.0001      0.000066    0.0001      0.000066    0.43     0.43     |
|                                                                                                              |
| RBEIMEvaluation                                                                                              |
|   rb_solve()                     1         0.0069      0.006882    0.0069      0.006882    22.67    22.67    |
|   read_offline_data_from_files() 1         0.0001      0.000092    0.0004      0.000438    0.30     1.44     |
|                                                                                                              |
| RBEvaluation                                                                                                 |
|   clear()                        3         0.0001      0.000020    0.0001      0.000020    0.19     0.19     |
|   compute_residual_dual_norm()   1         0.0005      0.000525    0.0005      0.000525    1.73     1.73     |
|   rb_solve()                     1         0.0001      0.000094    0.0075      0.007502    0.31     24.71    |
|   read_offline_data_from_files() 2         0.0007      0.000353    0.0008      0.000380    2.33     2.51     |
|   resize_data_structures()       2         0.0001      0.000027    0.0001      0.000027    0.18     0.18     |
 --------------------------------------------------------------------------------------------------------------
| Totals:                          3251      0.0304                                          100.00            |
 --------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running  ./reduced_basis_ex4-opt
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
