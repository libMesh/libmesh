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
  

<br><br>rbOOmit is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.
  

<br><br>You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


<br><br><h1>Reduced Basis Example 4 - Empirical Interpolation Method</h1>


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
        
          if(!online_mode)
          {
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
Read data from input file and print state
</div>

<div class ="fragment">
<pre>
            eim_construction.process_parameters_file(eim_parameters);
            eim_construction.print_info();
          
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
Set the rb_eval objects for the RBConstructions
</div>

<div class ="fragment">
<pre>
            eim_construction.rb_eval = &eim_rb_eval;
            rb_construction.rb_eval = &rb_eval;
          
        
</pre>
</div>
<div class = "comment">
Perform the EIM Greedy and write out the data
</div>

<div class ="fragment">
<pre>
            eim_construction.initialize_rb_construction();
            eim_construction.train_reduced_basis();
            eim_construction.rb_eval-&gt;write_offline_data_to_files("eim_data");
        
</pre>
</div>
<div class = "comment">
attach the EIM theta objects to the RBConstruction and RBEvaluation objects
</div>

<div class ="fragment">
<pre>
            eim_rb_eval.initialize_rb_theta_objects();
            rb_construction.rb_theta_expansion-&gt;attach_multiple_theta_q_f(eim_rb_eval.rb_eim_theta_vector);
            rb_construction.rb_eval-&gt;rb_theta_expansion-&gt;attach_multiple_theta_q_f(eim_rb_eval.rb_eim_theta_vector);
            
</pre>
</div>
<div class = "comment">
attach the EIM assembly objects to the RBConstruction object
</div>

<div class ="fragment">
<pre>
            eim_construction.initialize_EIM_F_objects();
            rb_construction.get_rb_assembly_expansion().attach_multiple_F_q_assembly(eim_construction.rb_eim_f_vector);
        
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
            rb_construction.rb_eval-&gt;write_offline_data_to_files("rb_data");
        
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
              eim_construction.rb_eval-&gt;write_out_basis_functions(eim_construction,"eim_data");
              rb_construction.rb_eval-&gt;write_out_basis_functions(rb_construction,"rb_data");
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
            eim_rb_eval.initialize_rb_theta_objects();
            rb_eval.rb_theta_expansion-&gt;attach_multiple_theta_q_f(eim_rb_eval.rb_eim_theta_vector);
            
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
compute the RB solution
</div>

<div class ="fragment">
<pre>
            rb_eval.set_current_parameters(online_mu_vector);
            rb_eval.print_current_parameters();
            rb_eval.rb_solve(online_N);
        
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
initialize the EquationSystems object by reading in the state that
was written out in the offline stage
</div>

<div class ="fragment">
<pre>
              equation_systems.read("equation_systems.dat", READ);
              RBConstruction& rb_construction =
                equation_systems.get_system&lt;RBConstruction&gt;("RB");
              RBConstruction& eim_construction =
                equation_systems.get_system&lt;RBConstruction&gt;("EIM");
              rb_construction.rb_eval = &rb_eval;
              eim_construction.rb_eval = &eim_rb_eval;
        
</pre>
</div>
<div class = "comment">
read in the data from files
</div>

<div class ="fragment">
<pre>
              eim_rb_eval.read_in_basis_functions(eim_construction,"eim_data");
              rb_eval.read_in_basis_functions(rb_construction,"rb_data");
        
</pre>
</div>
<div class = "comment">
load the EIM approximation for online_mu
</div>

<div class ="fragment">
<pre>
              eim_rb_eval.set_current_parameters(online_mu_vector);
              eim_rb_eval.rb_solve(eim_rb_eval.get_n_basis_functions());
        
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
    libmesh_example_assert(libMesh::default_solver_package() == PETSC_SOLVERS, <B><FONT COLOR="#BC8F8F">&quot;--enable-petsc&quot;</FONT></B>);
  
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
  
    SimpleRBEvaluation rb_eval;
  
    SimpleEIMEvaluation eim_rb_eval;
  
    <B><FONT COLOR="#A020F0">if</FONT></B>(!online_mode)
    {
      SimpleEIMConstruction &amp; eim_construction =
        equation_systems.add_system&lt;SimpleEIMConstruction&gt; (<B><FONT COLOR="#BC8F8F">&quot;EIM&quot;</FONT></B>);
      SimpleRBConstruction &amp; rb_construction =
        equation_systems.add_system&lt;SimpleRBConstruction&gt; (<B><FONT COLOR="#BC8F8F">&quot;RB&quot;</FONT></B>);
      
      equation_systems.init ();
    
      mesh.print_info();
      equation_systems.print_info();
  
      eim_construction.process_parameters_file(eim_parameters);
      eim_construction.print_info();
    
      rb_construction.process_parameters_file(rb_parameters);
      
      eim_construction.rb_eval = &amp;eim_rb_eval;
      rb_construction.rb_eval = &amp;rb_eval;
    
  
      eim_construction.initialize_rb_construction();
      eim_construction.train_reduced_basis();
      eim_construction.rb_eval-&gt;write_offline_data_to_files(<B><FONT COLOR="#BC8F8F">&quot;eim_data&quot;</FONT></B>);
  
      eim_rb_eval.initialize_rb_theta_objects();
      rb_construction.rb_theta_expansion-&gt;attach_multiple_theta_q_f(eim_rb_eval.rb_eim_theta_vector);
      rb_construction.rb_eval-&gt;rb_theta_expansion-&gt;attach_multiple_theta_q_f(eim_rb_eval.rb_eim_theta_vector);
      
      eim_construction.initialize_EIM_F_objects();
      rb_construction.get_rb_assembly_expansion().attach_multiple_F_q_assembly(eim_construction.rb_eim_f_vector);
  
      rb_construction.print_info();
  
      rb_construction.initialize_rb_construction();
      rb_construction.train_reduced_basis();
      rb_construction.rb_eval-&gt;write_offline_data_to_files(<B><FONT COLOR="#BC8F8F">&quot;rb_data&quot;</FONT></B>);
  
      <B><FONT COLOR="#A020F0">if</FONT></B>(store_basis_functions)
      {
        equation_systems.write(<B><FONT COLOR="#BC8F8F">&quot;equation_systems.dat&quot;</FONT></B>, WRITE);
  
        eim_construction.rb_eval-&gt;write_out_basis_functions(eim_construction,<B><FONT COLOR="#BC8F8F">&quot;eim_data&quot;</FONT></B>);
        rb_construction.rb_eval-&gt;write_out_basis_functions(rb_construction,<B><FONT COLOR="#BC8F8F">&quot;rb_data&quot;</FONT></B>);
      }
    }
    <B><FONT COLOR="#A020F0">else</FONT></B>
    {
      eim_rb_eval.read_offline_data_from_files(<B><FONT COLOR="#BC8F8F">&quot;eim_data&quot;</FONT></B>);
  
      eim_rb_eval.initialize_rb_theta_objects();
      rb_eval.rb_theta_expansion-&gt;attach_multiple_theta_q_f(eim_rb_eval.rb_eim_theta_vector);
      
      rb_eval.read_offline_data_from_files(<B><FONT COLOR="#BC8F8F">&quot;rb_data&quot;</FONT></B>);
  
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
  
      <B><FONT COLOR="#A020F0">if</FONT></B>(store_basis_functions)
      {
        equation_systems.read(<B><FONT COLOR="#BC8F8F">&quot;equation_systems.dat&quot;</FONT></B>, READ);
        RBConstruction&amp; rb_construction =
          equation_systems.get_system&lt;RBConstruction&gt;(<B><FONT COLOR="#BC8F8F">&quot;RB&quot;</FONT></B>);
        RBConstruction&amp; eim_construction =
          equation_systems.get_system&lt;RBConstruction&gt;(<B><FONT COLOR="#BC8F8F">&quot;EIM&quot;</FONT></B>);
        rb_construction.rb_eval = &amp;rb_eval;
        eim_construction.rb_eval = &amp;eim_rb_eval;
  
        eim_rb_eval.read_in_basis_functions(eim_construction,<B><FONT COLOR="#BC8F8F">&quot;eim_data&quot;</FONT></B>);
        rb_eval.read_in_basis_functions(rb_construction,<B><FONT COLOR="#BC8F8F">&quot;rb_data&quot;</FONT></B>);
  
        eim_rb_eval.set_current_parameters(online_mu_vector);
        eim_rb_eval.rb_solve(eim_rb_eval.get_n_basis_functions());
  
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
Compiling C++ (in optimized mode) reduced_basis_ex4.C...
Linking reduced_basis_ex4-opt...
***************************************************************
* Running  ./reduced_basis_ex4-opt
***************************************************************
 
*** Warning, This code is untested, experimental, or likely to see future API changes: src/reduced_basis/rb_parametrized_object.C, line 31, compiled Apr  7 2012 at 15:51:42 ***
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
Parameter 0: log scaling = 0
Parameter 1: log scaling = 0


RBConstruction parameters:
system name: EIM
constrained_problem: 0
Nmax: 20
Basis training error tolerance: 0.001
A_q operators attached: 0
F_q functions attached: 0
n_outputs: 0
Number of parameters: 2
Parameter 0: Min = -1, Max = 1
Parameter 1: Min = -1, Max = 1
n_training_samples: 25
single-matrix mode? 0
reuse preconditioner? 1
use a relative error bound in greedy? 0
write out data during basis training? 0
quiet mode? 1
parameter initialized to: 
mu[0] = -1
mu[1] = -1


RBEIMConstruction parameters:
best fit type: eim

Initializing training parameters with deterministic training set...
Parameter 0: log scaling = 0
Parameter 1: log scaling = 0


---- Performing Greedy basis enrichment ----

---- Basis dimension: 0 ----
Performing truth solve at parameter:
mu[0] = -1
mu[1] = -1

Enriching the RB space
Updating RB matrices

---- Basis dimension: 1 ----
Performing RB solves on training set
Maximum (absolute) error bound is 1.00428

Performing truth solve at parameter:
mu[0] = 1
mu[1] = 1

Enriching the RB space
Updating RB matrices

---- Basis dimension: 2 ----
Performing RB solves on training set
Maximum (absolute) error bound is 1.00428

Performing truth solve at parameter:
mu[0] = 1
mu[1] = -1

Enriching the RB space
Updating RB matrices

---- Basis dimension: 3 ----
Performing RB solves on training set
Maximum (absolute) error bound is 1.00428

Performing truth solve at parameter:
mu[0] = -1
mu[1] = 1

Enriching the RB space
Updating RB matrices

---- Basis dimension: 4 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.996306

Performing truth solve at parameter:
mu[0] = 0
mu[1] = 0

Enriching the RB space
Updating RB matrices

---- Basis dimension: 5 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.946406

Performing truth solve at parameter:
mu[0] = 0
mu[1] = -1

Enriching the RB space
Updating RB matrices

---- Basis dimension: 6 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.946117

Performing truth solve at parameter:
mu[0] = -1
mu[1] = 0

Enriching the RB space
Updating RB matrices

---- Basis dimension: 7 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.94264

Performing truth solve at parameter:
mu[0] = 0
mu[1] = 1

Enriching the RB space
Updating RB matrices

---- Basis dimension: 8 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.942187

Performing truth solve at parameter:
mu[0] = 1
mu[1] = 0

Enriching the RB space
Updating RB matrices

---- Basis dimension: 9 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.579353

Performing truth solve at parameter:
mu[0] = -0.5
mu[1] = -0.5

Enriching the RB space
Updating RB matrices

---- Basis dimension: 10 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.560329

Performing truth solve at parameter:
mu[0] = 0.5
mu[1] = -0.5

Enriching the RB space
Updating RB matrices

---- Basis dimension: 11 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.544525

Performing truth solve at parameter:
mu[0] = -0.5
mu[1] = 0.5

Enriching the RB space
Updating RB matrices

---- Basis dimension: 12 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.53201

Performing truth solve at parameter:
mu[0] = 0.5
mu[1] = 0.5

Enriching the RB space
Updating RB matrices

---- Basis dimension: 13 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.265185

Performing truth solve at parameter:
mu[0] = -0.5
mu[1] = -1

Enriching the RB space
Updating RB matrices

---- Basis dimension: 14 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.263071

Performing truth solve at parameter:
mu[0] = 0.5
mu[1] = 1

Enriching the RB space
Updating RB matrices

---- Basis dimension: 15 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.257056

Performing truth solve at parameter:
mu[0] = 1
mu[1] = -0.5

Enriching the RB space
Updating RB matrices

---- Basis dimension: 16 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.244274

Performing truth solve at parameter:
mu[0] = -1
mu[1] = 0.5

Enriching the RB space
Updating RB matrices

---- Basis dimension: 17 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.210069

Performing truth solve at parameter:
mu[0] = -1
mu[1] = -0.5

Enriching the RB space
Updating RB matrices

---- Basis dimension: 18 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.199985

Performing truth solve at parameter:
mu[0] = 0.5
mu[1] = -1

Enriching the RB space
Updating RB matrices

---- Basis dimension: 19 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.199126

Performing truth solve at parameter:
mu[0] = -0.5
mu[1] = 1

Enriching the RB space
Updating RB matrices

---- Basis dimension: 20 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.198364

Maximum number of basis functions reached: Nmax = 20.
Perform one more Greedy iteration for error bounds.
Performing truth solve at parameter:
mu[0] = 1
mu[1] = 0.5

Enriching the RB space
Updating RB matrices

---- Basis dimension: 20 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.198364

Extra Greedy iteration finished.

RBConstruction parameters:
system name: RB
constrained_problem: 0
Nmax: 15
Basis training error tolerance: 0.001
A_q operators attached: 1
F_q functions attached: 20
n_outputs: 0
Number of parameters: 2
Parameter 0: Min = -1, Max = 1
Parameter 1: Min = -1, Max = 1
n_training_samples: 100
single-matrix mode? 0
reuse preconditioner? 1
use a relative error bound in greedy? 0
write out data during basis training? 0
quiet mode? 1
parameter initialized to: 
mu[0] = -1
mu[1] = -1


---- Performing Greedy basis enrichment ----

---- Basis dimension: 0 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.387665

Performing truth solve at parameter:
mu[0] = 0.111111
mu[1] = 0.111111

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 1 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.158203

Performing truth solve at parameter:
mu[0] = 0.555556
mu[1] = -0.555556

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 2 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.155393

Performing truth solve at parameter:
mu[0] = -0.555556
mu[1] = -0.333333

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 3 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.126491

Performing truth solve at parameter:
mu[0] = -0.555556
mu[1] = 0.555556

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 4 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.0930572

Performing truth solve at parameter:
mu[0] = 0.555556
mu[1] = 0.777778

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 5 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.0562574

Performing truth solve at parameter:
mu[0] = 1
mu[1] = 0.111111

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 6 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.0465831

Performing truth solve at parameter:
mu[0] = -0.111111
mu[1] = -1

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 7 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.0339476

Performing truth solve at parameter:
mu[0] = -0.111111
mu[1] = 1

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 8 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.0236964

Performing truth solve at parameter:
mu[0] = -1
mu[1] = 1

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 9 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.0212792

Performing truth solve at parameter:
mu[0] = 1
mu[1] = -1

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 10 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.0201867

Performing truth solve at parameter:
mu[0] = -1
mu[1] = 0.111111

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 11 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.0134935

Performing truth solve at parameter:
mu[0] = 1
mu[1] = 1

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 12 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.0104616

Performing truth solve at parameter:
mu[0] = -1
mu[1] = -0.555556

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 13 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.0101337

Performing truth solve at parameter:
mu[0] = 0.333333
mu[1] = -1

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 14 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.00609637

Performing truth solve at parameter:
mu[0] = 1
mu[1] = -0.333333

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 15 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.00413525

Maximum number of basis functions reached: Nmax = 15
Writing out the basis functions...
Writing out the basis functions...

-------------------------------------------------------------------
| Time:           Sat Apr  7 16:02:22 2012                         |
| OS:             Linux                                            |
| HostName:       lkirk-home                                       |
| OS Release:     3.0.0-17-generic                                 |
| OS Version:     #30-Ubuntu SMP Thu Mar 8 20:45:39 UTC 2012       |
| Machine:        x86_64                                           |
| Username:       benkirk                                          |
| Configuration:  ./configure run on Sat Apr  7 15:49:27 CDT 2012  |
-------------------------------------------------------------------
 -----------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=27.0579, Active time=26.2887                                                    |
 -----------------------------------------------------------------------------------------------------------------
| Event                               nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                               w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|-----------------------------------------------------------------------------------------------------------------|
|                                                                                                                 |
|                                                                                                                 |
| DofMap                                                                                                          |
|   add_neighbors_to_send_list()      2         0.0003      0.000172    0.0003      0.000172    0.00     0.00     |
|   build_constraint_matrix()         13750     0.0113      0.000001    0.0113      0.000001    0.04     0.04     |
|   cnstrn_elem_mat_vec()             13750     0.0327      0.000002    0.0327      0.000002    0.12     0.12     |
|   compute_sparsity()                2         0.0030      0.001508    0.0040      0.001984    0.01     0.02     |
|   create_dof_constraints()          2         0.0015      0.000726    0.0019      0.000938    0.01     0.01     |
|   distribute_dofs()                 2         0.0006      0.000294    0.0022      0.001088    0.00     0.01     |
|   dof_indices()                     752315    0.4125      0.000001    0.4125      0.000001    1.57     1.57     |
|   prepare_send_list()               2         0.0000      0.000000    0.0000      0.000000    0.00     0.00     |
|   reinit()                          2         0.0016      0.000792    0.0016      0.000792    0.01     0.01     |
|                                                                                                                 |
| EquationSystems                                                                                                 |
|   write()                           1         0.0183      0.018339    0.0187      0.018696    0.07     0.07     |
|                                                                                                                 |
| FE                                                                                                              |
|   compute_affine_map()              371050    0.2945      0.000001    0.2945      0.000001    1.12     1.12     |
|   compute_face_map()                2300      0.0068      0.000003    0.0137      0.000006    0.03     0.05     |
|   compute_shape_functions()         371050    0.1141      0.000000    0.1141      0.000000    0.43     0.43     |
|   init_face_shape_functions()       23        0.0001      0.000002    0.0001      0.000002    0.00     0.00     |
|   init_shape_functions()            2890      0.0148      0.000005    0.0148      0.000005    0.06     0.06     |
|   inverse_map()                     55480     0.0726      0.000001    0.0726      0.000001    0.28     0.28     |
|                                                                                                                 |
| Mesh                                                                                                            |
|   find_neighbors()                  1         0.0011      0.001112    0.0011      0.001112    0.00     0.00     |
|   renumber_nodes_and_elem()         2         0.0001      0.000037    0.0001      0.000037    0.00     0.00     |
|                                                                                                                 |
| MeshCommunication                                                                                               |
|   assign_global_indices()           1         0.0151      0.015055    0.0151      0.015066    0.06     0.06     |
|                                                                                                                 |
| MeshTools::Generation                                                                                           |
|   build_cube()                      1         0.0007      0.000716    0.0007      0.000716    0.00     0.00     |
|                                                                                                                 |
| Parallel                                                                                                        |
|   allgather()                       6         0.0000      0.000001    0.0000      0.000001    0.00     0.00     |
|   receive()                         148       0.0004      0.000003    0.0004      0.000003    0.00     0.00     |
|   send()                            148       0.0009      0.000006    0.0009      0.000006    0.00     0.00     |
|   send_receive()                    4         0.0000      0.000002    0.0000      0.000002    0.00     0.00     |
|                                                                                                                 |
| Partitioner                                                                                                     |
|   single_partition()                1         0.0001      0.000072    0.0001      0.000072    0.00     0.00     |
|                                                                                                                 |
| PetscLinearSolver                                                                                               |
|   solve()                           596       0.3744      0.000628    0.3744      0.000628    1.42     1.42     |
|                                                                                                                 |
| PointLocatorTree                                                                                                |
|   init(no master)                   1         0.0010      0.000965    0.0010      0.000965    0.00     0.00     |
|   operator()                        440       0.0050      0.000011    0.0059      0.000013    0.02     0.02     |
|                                                                                                                 |
| RBConstruction                                                                                                  |
|   add_scaled_matrix_and_vector()    23        0.1179      0.005125    0.4281      0.018611    0.45     1.63     |
|   clear()                           2         0.0005      0.000245    0.0005      0.000245    0.00     0.00     |
|   compute_Fq_representor_norms()    2         0.0097      0.004836    0.0586      0.029291    0.04     0.22     |
|   compute_max_error_bound()         37        0.0068      0.000185    25.7870     0.696946    0.03     98.09    |
|   enrich_RB_space()                 15        0.0061      0.000404    0.0061      0.000404    0.02     0.02     |
|   train_reduced_basis()             2         0.0045      0.002266    26.2311     13.115566   0.02     99.78    |
|   truth_assembly()                  15        0.0109      0.000725    0.0189      0.001261    0.04     0.07     |
|   truth_solve()                     15        0.0010      0.000068    0.0650      0.004332    0.00     0.25     |
|   update_RB_system_matrices()       36        0.0173      0.000479    0.0173      0.000479    0.07     0.07     |
|   update_residual_terms()           15        0.0252      0.001682    0.0715      0.004767    0.10     0.27     |
|                                                                                                                 |
| RBEIMConstruction                                                                                               |
|   compute_best_fit_error()          525       0.0158      0.000030    3.0142      0.005741    0.06     11.47    |
|   enrich_RB_space()                 21        0.0505      0.002404    0.0930      0.004429    0.19     0.35     |
|   evaluate_current_basis_function() 12500     0.0919      0.000007    0.1845      0.000015    0.35     0.70     |
|   set_current_basis_function()      12500     0.0044      0.000000    0.0044      0.000000    0.02     0.02     |
|   truth_solve()                     546       1.9062      0.003491    3.1077      0.005692    7.25     11.82    |
|   update_RB_system_matrices()       21        0.0042      0.000201    0.0150      0.000712    0.02     0.06     |
|                                                                                                                 |
| RBEIMEvaluation                                                                                                 |
|   rb_solve()                        944845    21.7485     0.000023    21.7485     0.000023    82.73    82.73    |
|   write_offline_data_to_files()     1         0.0006      0.000561    0.0015      0.001477    0.00     0.01     |
|                                                                                                                 |
| RBEvaluation                                                                                                    |
|   clear()                           3         0.0001      0.000034    0.0001      0.000044    0.00     0.00     |
|   clear_riesz_representors()        5         0.0000      0.000006    0.0000      0.000006    0.00     0.00     |
|   compute_residual_dual_norm()      1600      0.8307      0.000519    21.9691     0.013731    3.16     83.57    |
|   rb_solve()                        1600      0.0487      0.000030    22.7656     0.014229    0.19     86.60    |
|   resize_data_structures()          2         0.0001      0.000033    0.0001      0.000034    0.00     0.00     |
|   write_offline_data_to_files()     2         0.0037      0.001873    0.0037      0.001873    0.01     0.01     |
 -----------------------------------------------------------------------------------------------------------------
| Totals:                             2558303   26.2887                                         100.00            |
 -----------------------------------------------------------------------------------------------------------------

*** Warning, This code is untested, experimental, or likely to see future API changes: src/reduced_basis/rb_parametrized_object.C, line 31, compiled Apr  7 2012 at 15:51:42 ***
mu[0] = -0.6
mu[1] = 0.7

Reading in the basis functions...
Finished reading in the basis functions...
Reading in the basis functions...
Finished reading in the basis functions...

-------------------------------------------------------------------
| Time:           Sat Apr  7 16:02:23 2012                         |
| OS:             Linux                                            |
| HostName:       lkirk-home                                       |
| OS Release:     3.0.0-17-generic                                 |
| OS Version:     #30-Ubuntu SMP Thu Mar 8 20:45:39 UTC 2012       |
| Machine:        x86_64                                           |
| Username:       benkirk                                          |
| Configuration:  ./configure run on Sat Apr  7 15:49:27 CDT 2012  |
-------------------------------------------------------------------
 --------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=0.382762, Active time=0.076274                                               |
 --------------------------------------------------------------------------------------------------------------
| Event                            nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                            w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|--------------------------------------------------------------------------------------------------------------|
|                                                                                                              |
|                                                                                                              |
| DofMap                                                                                                       |
|   add_neighbors_to_send_list()   2         0.0005      0.000228    0.0005      0.000228    0.60     0.60     |
|   compute_sparsity()             2         0.0039      0.001929    0.0051      0.002534    5.06     6.64     |
|   create_dof_constraints()       2         0.0004      0.000218    0.0004      0.000218    0.57     0.57     |
|   distribute_dofs()              2         0.0007      0.000363    0.0028      0.001418    0.95     3.72     |
|   dof_indices()                  2500      0.0017      0.000001    0.0017      0.000001    2.17     2.17     |
|   prepare_send_list()            2         0.0000      0.000002    0.0000      0.000002    0.00     0.00     |
|   reinit()                       2         0.0021      0.001053    0.0021      0.001053    2.76     2.76     |
|                                                                                                              |
| EquationSystems                                                                                              |
|   build_solution_vector()        1         0.0014      0.001419    0.0024      0.002422    1.86     3.18     |
|   read()                         1         0.0192      0.019223    0.0420      0.042042    25.20    55.12    |
|   update()                       1         0.0001      0.000091    0.0001      0.000091    0.12     0.12     |
|                                                                                                              |
| ExodusII_IO                                                                                                  |
|   write_nodal_data()             1         0.0017      0.001726    0.0017      0.001726    2.26     2.26     |
|                                                                                                              |
| Mesh                                                                                                         |
|   find_neighbors()               1         0.0011      0.001105    0.0011      0.001105    1.45     1.45     |
|   renumber_nodes_and_elem()      2         0.0001      0.000038    0.0001      0.000038    0.10     0.10     |
|                                                                                                              |
| MeshCommunication                                                                                            |
|   assign_global_indices()        1         0.0140      0.013958    0.0140      0.013968    18.30    18.31    |
|                                                                                                              |
| MeshOutput                                                                                                   |
|   write_equation_systems()       1         0.0000      0.000036    0.0042      0.004184    0.05     5.49     |
|                                                                                                              |
| MeshTools::Generation                                                                                        |
|   build_cube()                   1         0.0006      0.000580    0.0006      0.000580    0.76     0.76     |
|                                                                                                              |
| Parallel                                                                                                     |
|   allgather()                    80        0.0001      0.000001    0.0001      0.000001    0.13     0.13     |
|   send_receive()                 4         0.0000      0.000001    0.0000      0.000001    0.01     0.01     |
|                                                                                                              |
| Partitioner                                                                                                  |
|   single_partition()             1         0.0001      0.000070    0.0001      0.000070    0.09     0.09     |
|                                                                                                              |
| RBConstruction                                                                                               |
|   clear()                        4         0.0002      0.000053    0.0002      0.000053    0.28     0.28     |
|   load_rb_solution()             2         0.0002      0.000079    0.0002      0.000079    0.21     0.21     |
|                                                                                                              |
| RBEIMEvaluation                                                                                              |
|   rb_solve()                     741       0.0236      0.000032    0.0236      0.000032    30.93    30.93    |
|   read_offline_data_from_files() 1         0.0006      0.000561    0.0013      0.001290    0.74     1.69     |
|                                                                                                              |
| RBEvaluation                                                                                                 |
|   clear()                        3         0.0001      0.000030    0.0001      0.000031    0.12     0.12     |
|   clear_riesz_representors()     5         0.0000      0.000001    0.0000      0.000001    0.01     0.01     |
|   compute_residual_dual_norm()   1         0.0010      0.000968    0.0236      0.023584    1.27     30.92    |
|   rb_solve()                     1         0.0001      0.000090    0.0247      0.024727    0.12     32.42    |
|   read_offline_data_from_files() 2         0.0029      0.001433    0.0030      0.001488    3.76     3.90     |
|   resize_data_structures()       2         0.0001      0.000051    0.0001      0.000054    0.14     0.14     |
 --------------------------------------------------------------------------------------------------------------
| Totals:                          3369      0.0763                                          100.00            |
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
