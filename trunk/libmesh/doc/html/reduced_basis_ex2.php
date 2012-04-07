<?php $root=""; ?>
<?php require($root."navigation.php"); ?>
<html>
<head>
  <?php load_style($root); ?>
</head>
 
<body>
 
<?php make_navigation("reduced_basis_ex2",$root)?>
 
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
<h1>Reduced Basis Example 2 - Successive Constraint Method</h1>


<br><br>In this example we extend reduced_basis_ex1 to solve a steady convection-diffusion
problem on the unit square via the Reduced Basis Method. In this case, we modify the
PDE so that it no longer has a parameter-independent coercivity constant. Therefore,
in order to obtain an error bound, we need to employ the Successive Constraint
Method (SCM) implemented in RBSCMConstruction/RBSCMEvaluation to obtain a
parameter-dependent lower bound for the coercivity constant.


<br><br>The PDE being solved is div(k*grad(u)) + Beta*grad(u) = f
k is the diffusion coefficient :
- constant in the domain 0<=x<0.5 , its value is given by the first parameter mu[0]
- constant in the domain 0.5<=x<=1 , its value is given by the second parameter mu[1]
Beta is the convection velocity :
- constant in the whole domain
- equal to zero in the y-direction
- its value in the x-direction is given by the third (and last) parameter mu[2]
Boundary conditions :
- dyu=0 on top and bottom
- u=0 on the left side
- dxu + Beta*u = 0 on the right side


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
This example requires SLEPc and GLPK
</div>

<div class ="fragment">
<pre>
        #if !defined(LIBMESH_HAVE_SLEPC) || !defined(LIBMESH_HAVE_GLPK)
          libmesh_example_assert(false, "--enable-slepc --enable-glpk");
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
Parse the input file (reduced_basis_ex2.in) using GetPot
</div>

<div class ="fragment">
<pre>
          std::string parameters_filename = "reduced_basis_ex2.in";
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
        
</pre>
</div>
<div class = "comment">
We also need a SCM evaluation object to perform SCM calculations
</div>

<div class ="fragment">
<pre>
          RBSCMEvaluation rb_scm_eval;
          rb_scm_eval.rb_theta_expansion = rb_eval.rb_theta_expansion;
        
</pre>
</div>
<div class = "comment">
Tell rb_eval about rb_scm_eval
</div>

<div class ="fragment">
<pre>
          rb_eval.rb_scm_eval = &rb_scm_eval;
          
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
Initialize the SCM Construction object
</div>

<div class ="fragment">
<pre>
            RBSCMConstruction & rb_scm_con =
              equation_systems.add_system&lt;RBSCMConstruction&gt; ("RBSCMConvectionDiffusion");
            rb_scm_con.set_RB_system_name("RBConvectionDiffusion");
            rb_scm_con.add_variable("p", FIRST);
        
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
Set parameters for the eigenvalue problems that will be solved by rb_scm_con
</div>

<div class ="fragment">
<pre>
            equation_systems.parameters.set&lt;unsigned int&gt;("eigenpairs")    = 1;
            equation_systems.parameters.set&lt;unsigned int&gt;("basis vectors") = 3;
            equation_systems.parameters.set&lt;unsigned int&gt;
              ("linear solver maximum iterations") = 1000;
        
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
            rb_scm_con.process_parameters_file(parameters_filename);
        
</pre>
</div>
<div class = "comment">
Need to give rb_scm_con a pointer to the theta expansion
</div>

<div class ="fragment">
<pre>
            rb_scm_con.rb_theta_expansion = rb_con.rb_theta_expansion;
          
</pre>
</div>
<div class = "comment">
Finally, need to give rb_scm_con and rb_eval a pointer to the
SCM evaluation object, rb_scm_eval
</div>

<div class ="fragment">
<pre>
            rb_scm_con.rb_scm_eval = &rb_scm_eval;
        
</pre>
</div>
<div class = "comment">
Print out info that describes the current setup of rb_con
</div>

<div class ="fragment">
<pre>
            rb_con.print_info();
            rb_scm_con.print_info();
        
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
Perform the SCM Greedy algorithm to derive the data required
for rb_scm_eval to provide a coercivity lower bound.
</div>

<div class ="fragment">
<pre>
            rb_scm_con.perform_SCM_greedy();
        
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
            rb_scm_con.rb_scm_eval-&gt;write_offline_data_to_files();
            
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
            rb_scm_eval.set_current_parameters(online_mu_vector);
            rb_eval.print_current_parameters();
            
</pre>
</div>
<div class = "comment">
Read in the reduced basis data
</div>

<div class ="fragment">
<pre>
            rb_eval.read_offline_data_from_files();
            rb_scm_eval.read_offline_data_from_files();   
         
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
        
        #endif // LIBMESH_HAVE_SLEPC && LIBMESH_HAVE_GLPK
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
  
  #<B><FONT COLOR="#A020F0">if</FONT></B> !defined(LIBMESH_HAVE_SLEPC) || !defined(LIBMESH_HAVE_GLPK)
    libmesh_example_assert(false, <B><FONT COLOR="#BC8F8F">&quot;--enable-slepc --enable-glpk&quot;</FONT></B>);
  #<B><FONT COLOR="#A020F0">else</FONT></B>
  
  #<B><FONT COLOR="#A020F0">if</FONT></B> !defined(LIBMESH_HAVE_XDR)
    libmesh_example_assert(false, <B><FONT COLOR="#BC8F8F">&quot;--enable-xdr&quot;</FONT></B>);
  #elif defined(LIBMESH_DEFAULT_SINGLE_PRECISION)
    libmesh_example_assert(false, <B><FONT COLOR="#BC8F8F">&quot;--disable-singleprecision&quot;</FONT></B>);
  #endif
    libmesh_example_assert(libMesh::default_solver_package() == PETSC_SOLVERS, <B><FONT COLOR="#BC8F8F">&quot;--enable-petsc&quot;</FONT></B>);
  
    libmesh_example_assert(2 &lt;= LIBMESH_DIM, <B><FONT COLOR="#BC8F8F">&quot;2D support&quot;</FONT></B>);
    
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::string parameters_filename = <B><FONT COLOR="#BC8F8F">&quot;reduced_basis_ex2.in&quot;</FONT></B>;
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
  
    RBSCMEvaluation rb_scm_eval;
    rb_scm_eval.rb_theta_expansion = rb_eval.rb_theta_expansion;
  
    rb_eval.rb_scm_eval = &amp;rb_scm_eval;
    
    <B><FONT COLOR="#A020F0">if</FONT></B>(!online_mode) <I><FONT COLOR="#B22222">// Perform the Offline stage of the RB method
</FONT></I>    {
      SimpleRBConstruction &amp; rb_con =
        equation_systems.add_system&lt;SimpleRBConstruction&gt; (<B><FONT COLOR="#BC8F8F">&quot;RBConvectionDiffusion&quot;</FONT></B>);
  
      RBSCMConstruction &amp; rb_scm_con =
        equation_systems.add_system&lt;RBSCMConstruction&gt; (<B><FONT COLOR="#BC8F8F">&quot;RBSCMConvectionDiffusion&quot;</FONT></B>);
      rb_scm_con.set_RB_system_name(<B><FONT COLOR="#BC8F8F">&quot;RBConvectionDiffusion&quot;</FONT></B>);
      rb_scm_con.add_variable(<B><FONT COLOR="#BC8F8F">&quot;p&quot;</FONT></B>, FIRST);
  
      equation_systems.init ();
  
      equation_systems.print_info();
      mesh.print_info();
  
      equation_systems.parameters.set&lt;<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B>&gt;(<B><FONT COLOR="#BC8F8F">&quot;eigenpairs&quot;</FONT></B>)    = 1;
      equation_systems.parameters.set&lt;<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B>&gt;(<B><FONT COLOR="#BC8F8F">&quot;basis vectors&quot;</FONT></B>) = 3;
      equation_systems.parameters.set&lt;<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B>&gt;
        (<B><FONT COLOR="#BC8F8F">&quot;linear solver maximum iterations&quot;</FONT></B>) = 1000;
  
      rb_con.rb_eval = &amp;rb_eval;
  
      rb_con.process_parameters_file(parameters_filename);
      rb_scm_con.process_parameters_file(parameters_filename);
  
      rb_scm_con.rb_theta_expansion = rb_con.rb_theta_expansion;
    
      rb_scm_con.rb_scm_eval = &amp;rb_scm_eval;
  
      rb_con.print_info();
      rb_scm_con.print_info();
  
      rb_con.initialize_rb_construction();
      
      rb_scm_con.perform_SCM_greedy();
  
      rb_con.train_reduced_basis();
      
      rb_con.rb_eval-&gt;write_offline_data_to_files();
      rb_scm_con.rb_scm_eval-&gt;write_offline_data_to_files();
      
      <B><FONT COLOR="#A020F0">if</FONT></B>(store_basis_functions)
      {
        equation_systems.write(<B><FONT COLOR="#BC8F8F">&quot;equation_systems.dat&quot;</FONT></B>, WRITE);
        
        rb_con.rb_eval-&gt;write_out_basis_functions(rb_con);
      }
    }
    <B><FONT COLOR="#A020F0">else</FONT></B> <I><FONT COLOR="#B22222">// Perform the Online stage of the RB method
</FONT></I>    {
  
      <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> online_N = infile(<B><FONT COLOR="#BC8F8F">&quot;online_N&quot;</FONT></B>,1);
      <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_parameters = infile(<B><FONT COLOR="#BC8F8F">&quot;n_parameters&quot;</FONT></B>,1);
      <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;Real&gt; online_mu_vector(n_parameters);
      <B><FONT COLOR="#A020F0">for</FONT></B>(<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i&lt;n_parameters; i++)
      {
        online_mu_vector[i] = infile(<B><FONT COLOR="#BC8F8F">&quot;online_mu&quot;</FONT></B>, online_mu_vector[i], i);
      }
  
      rb_eval.set_current_parameters(online_mu_vector);
      rb_scm_eval.set_current_parameters(online_mu_vector);
      rb_eval.print_current_parameters();
      
      rb_eval.read_offline_data_from_files();
      rb_scm_eval.read_offline_data_from_files();   
   
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
  
  #endif <I><FONT COLOR="#B22222">// LIBMESH_HAVE_SLEPC &amp;&amp; LIBMESH_HAVE_GLPK
</FONT></I>  }
  
</pre> 
<a name="output"></a> 
<br><br><br> <h1> The console output of the program: </h1> 
<pre>
Updated .depend
Compiling C++ (in optimized mode) reduced_basis_ex2.C...
Linking reduced_basis_ex2-opt...
***************************************************************
* Running  ./reduced_basis_ex2-opt
***************************************************************
 
*** Warning, This code is untested, experimental, or likely to see future API changes: src/reduced_basis/rb_parametrized_object.C, line 31, compiled Apr  7 2012 at 15:51:42 ***
 EquationSystems
  n_systems()=2
   System #0, "RBConvectionDiffusion"
    Type "RBConstruction"
    Variables="u" 
    Finite Element Types="LAGRANGE", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" 
    Approximation Orders="FIRST", "THIRD" 
    n_dofs()=676
    n_local_dofs()=676
    n_constrained_dofs()=26
    n_local_constrained_dofs()=26
    n_vectors()=1
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 8.54438
      Average Off-Processor Bandwidth <= 0
      Maximum  On-Processor Bandwidth <= 9
      Maximum Off-Processor Bandwidth <= 0
    DofMap Constraints
      Number of DoF Constraints = 26
      Average DoF Constraint Length= 0
      Number of Node Constraints = 0
   System #1, "RBSCMConvectionDiffusion"
    Type "Eigen"
    Variables="p" 
    Finite Element Types="LAGRANGE", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" 
    Approximation Orders="FIRST", "THIRD" 
    n_dofs()=676
    n_local_dofs()=676
    n_constrained_dofs()=0
    n_local_constrained_dofs()=0
    n_vectors()=0
    n_matrices()=2
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 8.54438
      Average Off-Processor Bandwidth <= 0
      Maximum  On-Processor Bandwidth <= 9
      Maximum Off-Processor Bandwidth <= 0
    DofMap Constraints
      Number of DoF Constraints = 0
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

Initializing training parameters with random training set...
Parameter 0: log scaling = 1
Parameter 1: log scaling = 1
Parameter 2: log scaling = 1

Initializing training parameters with random training set...
Parameter 0: log scaling = 0
Parameter 1: log scaling = 0
Parameter 2: log scaling = 0


RBConstruction parameters:
system name: RBConvectionDiffusion
constrained_problem: 0
Nmax: 10
Basis training error tolerance: 1e-05
A_q operators attached: 3
F_q functions attached: 1
n_outputs: 4
output 0, Q_l = 1
output 1, Q_l = 1
output 2, Q_l = 1
output 3, Q_l = 1
Number of parameters: 3
Parameter 0: Min = 0.1, Max = 1
Parameter 1: Min = 0.1, Max = 1
Parameter 2: Min = 0.01, Max = 0.1
n_training_samples: 100
single-matrix mode? 0
reuse preconditioner? 1
use a relative error bound in greedy? 0
write out data during basis training? 0
quiet mode? 1
parameter initialized to: 
mu[0] = 0.1
mu[1] = 0.1
mu[2] = 0.01


RBSCMConstruction parameters:
system name: RBSCMConvectionDiffusion
SCM Greedy tolerance: 0.1
A_q operators attached: 3
Number of parameters: 3
Parameter 0: Min = 0.1, Max = 1
Parameter 1: Min = 0.1, Max = 1
Parameter 2: Min = 0.01, Max = 0.1
n_training_samples: 100


B_min(0) = -1.60614e-15
B_max(0) = 0.999932

B_min(1) = -1.10495e-15
B_max(1) = 0.999933

B_min(2) = -0.380786
B_max(2) = 1.0479e-16

SCM: Added mu = (0.856169,0.534242,0.0981491)

Stability constant for C_J(0) = 0.434784

SCM iteration 0, max_SCM_error = 0.877807

SCM: Added mu = (0.156786,0.83771,0.0281991)

-----------------------------------


Stability constant for C_J(1) = 0.112264

SCM iteration 1, max_SCM_error = 0.407185

SCM: Added mu = (0.157754,0.162916,0.0757656)

-----------------------------------


Stability constant for C_J(2) = 0.0890664

SCM iteration 2, max_SCM_error = 0.369496

SCM: Added mu = (0.911987,0.146745,0.0687817)

-----------------------------------


Stability constant for C_J(3) = 0.116078

SCM iteration 3, max_SCM_error = 0.228715

SCM: Added mu = (0.83329,0.932839,0.0283193)

-----------------------------------


Stability constant for C_J(4) = 0.595895

SCM iteration 4, max_SCM_error = 0.125248

SCM: Added mu = (0.927124,0.759389,0.0136778)

-----------------------------------


Stability constant for C_J(5) = 0.610675

SCM iteration 5, max_SCM_error = 0.103702

SCM: Added mu = (0.495604,0.131879,0.0948746)

-----------------------------------


Stability constant for C_J(6) = 0.0944009

SCM iteration 6, max_SCM_error = 0.0838038

SCM tolerance of 0.1 reached.

Compute output dual norms
output_dual_norms[0][0] = 0.839698
output_dual_norms[1][0] = 0.318298
output_dual_norms[2][0] = 0.318298
output_dual_norms[3][0] = 0.839698

---- Performing Greedy basis enrichment ----

---- Basis dimension: 0 ----
Performing RB solves on training set
Maximum (absolute) error bound is 7.47733

Performing truth solve at parameter:
mu[0] = 0.115923
mu[1] = 0.117464
mu[2] = 0.0537934

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 1 ----
Performing RB solves on training set
Maximum (absolute) error bound is 7.41026

Performing truth solve at parameter:
mu[0] = 0.115637
mu[1] = 0.660203
mu[2] = 0.0159299

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 2 ----
Performing RB solves on training set
Maximum (absolute) error bound is 1.98631

Performing truth solve at parameter:
mu[0] = 0.798377
mu[1] = 0.112704
mu[2] = 0.0449915

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 3 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.812427

Performing truth solve at parameter:
mu[0] = 0.275144
mu[1] = 0.108498
mu[2] = 0.0877104

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 4 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.0313249

Performing truth solve at parameter:
mu[0] = 0.121915
mu[1] = 0.575089
mu[2] = 0.043415

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 5 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.0153797

Performing truth solve at parameter:
mu[0] = 0.46529
mu[1] = 0.113255
mu[2] = 0.0124091

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 6 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.00132377

Performing truth solve at parameter:
mu[0] = 0.169421
mu[1] = 0.339475
mu[2] = 0.0909153

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 7 ----
Performing RB solves on training set
Maximum (absolute) error bound is 2.62726e-05

Performing truth solve at parameter:
mu[0] = 0.109466
mu[1] = 0.118121
mu[2] = 0.0229933

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 8 ----
Performing RB solves on training set
Maximum (absolute) error bound is 3.84317e-06

Specified error tolerance reached.
In RBSCMEvaluation::write_offline_data_to_files, directory offline_data already exists, overwriting contents.
Writing out the basis functions...

-------------------------------------------------------------------
| Time:           Sat Apr  7 16:01:24 2012                         |
| OS:             Linux                                            |
| HostName:       lkirk-home                                       |
| OS Release:     3.0.0-17-generic                                 |
| OS Version:     #30-Ubuntu SMP Thu Mar 8 20:45:39 UTC 2012       |
| Machine:        x86_64                                           |
| Username:       benkirk                                          |
| Configuration:  ./configure run on Sat Apr  7 15:49:27 CDT 2012  |
-------------------------------------------------------------------
 --------------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=13.6448, Active time=13.4575                                                       |
 --------------------------------------------------------------------------------------------------------------------
| Event                                  nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                                  w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|--------------------------------------------------------------------------------------------------------------------|
|                                                                                                                    |
|                                                                                                                    |
| CondensedEigenSystem                                                                                               |
|   get_eigenpair()                      13        0.0029      0.000220    0.0029      0.000220    0.02     0.02     |
|   solve()                              13        0.0090      0.000689    12.5368     0.964366    0.07     93.16    |
|                                                                                                                    |
| DofMap                                                                                                             |
|   add_neighbors_to_send_list()         2         0.0004      0.000178    0.0004      0.000178    0.00     0.00     |
|   build_constraint_matrix()            33750     0.0205      0.000001    0.0205      0.000001    0.15     0.15     |
|   cnstrn_elem_mat_vec()                33750     0.0298      0.000001    0.0298      0.000001    0.22     0.22     |
|   compute_sparsity()                   2         0.0028      0.001379    0.0036      0.001806    0.02     0.03     |
|   create_dof_constraints()             2         0.0014      0.000709    0.0019      0.000928    0.01     0.01     |
|   distribute_dofs()                    2         0.0006      0.000287    0.0021      0.001072    0.00     0.02     |
|   dof_indices()                        69375     0.0585      0.000001    0.0585      0.000001    0.43     0.43     |
|   prepare_send_list()                  2         0.0000      0.000001    0.0000      0.000001    0.00     0.00     |
|   reinit()                             2         0.0016      0.000783    0.0016      0.000783    0.01     0.01     |
|                                                                                                                    |
| EquationSystems                                                                                                    |
|   write()                              1         0.0197      0.019688    0.0199      0.019853    0.15     0.15     |
|                                                                                                                    |
| FE                                                                                                                 |
|   compute_affine_map()                 39150     0.0367      0.000001    0.0367      0.000001    0.27     0.27     |
|   compute_face_map()                   5400      0.0201      0.000004    0.0418      0.000008    0.15     0.31     |
|   compute_shape_functions()            39150     0.0250      0.000001    0.0250      0.000001    0.19     0.19     |
|   init_face_shape_functions()          54        0.0002      0.000003    0.0002      0.000003    0.00     0.00     |
|   init_shape_functions()               5454      0.0267      0.000005    0.0267      0.000005    0.20     0.20     |
|   inverse_map()                        10800     0.0195      0.000002    0.0195      0.000002    0.14     0.14     |
|                                                                                                                    |
| Mesh                                                                                                               |
|   find_neighbors()                     1         0.0012      0.001151    0.0012      0.001151    0.01     0.01     |
|   renumber_nodes_and_elem()            2         0.0001      0.000035    0.0001      0.000035    0.00     0.00     |
|                                                                                                                    |
| MeshCommunication                                                                                                  |
|   assign_global_indices()              1         0.0142      0.014212    0.0142      0.014224    0.11     0.11     |
|                                                                                                                    |
| MeshTools::Generation                                                                                              |
|   build_cube()                         1         0.0010      0.001006    0.0010      0.001006    0.01     0.01     |
|                                                                                                                    |
| Parallel                                                                                                           |
|   allgather()                          6         0.0000      0.000000    0.0000      0.000000    0.00     0.00     |
|   receive()                            40        0.0001      0.000004    0.0001      0.000004    0.00     0.00     |
|   send()                               40        0.0003      0.000007    0.0003      0.000007    0.00     0.00     |
|   send_receive()                       4         0.0000      0.000002    0.0000      0.000002    0.00     0.00     |
|                                                                                                                    |
| Partitioner                                                                                                        |
|   single_partition()                   1         0.0001      0.000069    0.0001      0.000069    0.00     0.00     |
|                                                                                                                    |
| PetscLinearSolver                                                                                                  |
|   solve()                              37        0.1584      0.004281    0.1584      0.004281    1.18     1.18     |
|                                                                                                                    |
| RBConstruction                                                                                                     |
|   add_scaled_Aq()                      45        0.0012      0.000027    0.5473      0.012163    0.01     4.07     |
|   add_scaled_matrix_and_vector()       54        0.3431      0.006354    0.6250      0.011574    2.55     4.64     |
|   clear()                              1         0.0004      0.000386    0.0004      0.000386    0.00     0.00     |
|   compute_Fq_representor_norms()       1         0.0004      0.000358    0.0043      0.004325    0.00     0.03     |
|   compute_max_error_bound()            9         0.0024      0.000270    0.0484      0.005377    0.02     0.36     |
|   compute_output_dual_norms()          1         0.0111      0.011052    0.0297      0.029657    0.08     0.22     |
|   enrich_RB_space()                    8         0.0019      0.000242    0.0019      0.000242    0.01     0.01     |
|   train_reduced_basis()                1         0.0011      0.001144    0.2530      0.253016    0.01     1.88     |
|   truth_assembly()                     8         0.0069      0.000856    0.0069      0.000856    0.05     0.05     |
|   truth_solve()                        8         0.0007      0.000085    0.0441      0.005513    0.01     0.33     |
|   update_RB_system_matrices()          8         0.0066      0.000830    0.0066      0.000830    0.05     0.05     |
|   update_residual_terms()              8         0.0175      0.002184    0.1167      0.014593    0.13     0.87     |
|                                                                                                                    |
| RBEvaluation                                                                                                       |
|   clear()                              1         0.0001      0.000058    0.0001      0.000117    0.00     0.00     |
|   clear_riesz_representors()           2         0.0001      0.000029    0.0001      0.000029    0.00     0.00     |
|   compute_residual_dual_norm()         900       0.0053      0.000006    0.0053      0.000006    0.04     0.04     |
|   rb_solve()                           900       0.0086      0.000010    0.0458      0.000051    0.06     0.34     |
|   resize_data_structures()             1         0.0001      0.000068    0.0001      0.000069    0.00     0.00     |
|   write_offline_data_to_files()        1         0.0016      0.001591    0.0016      0.001591    0.01     0.01     |
|                                                                                                                    |
| RBSCMConstruction                                                                                                  |
|   add_scaled_symm_Aq()                 45        0.0002      0.000004    0.5475      0.012167    0.00     4.07     |
|   compute_SCM_bounding_box()           1         0.0009      0.000909    5.3935      5.393472    0.01     40.08    |
|   compute_SCM_bounds_on_training_set() 7         0.0028      0.000401    0.0246      0.003519    0.02     0.18     |
|   enrich_C_J()                         7         0.0002      0.000028    0.0002      0.000028    0.00     0.00     |
|   evaluate_stability_constant()        7         0.0037      0.000527    7.6983      1.099757    0.03     57.20    |
|   perform_SCM_greedy()                 1         0.0088      0.008844    13.1255     13.125455   0.07     97.53    |
|                                                                                                                    |
| RBSCMEvaluation                                                                                                    |
|   get_SCM_LB()                         1600      0.0525      0.000033    0.0525      0.000033    0.39     0.39     |
|   get_SCM_UB()                         700       0.0008      0.000001    0.0008      0.000001    0.01     0.01     |
|   write_offline_data_to_files()        1         0.0002      0.000230    0.0002      0.000230    0.00     0.00     |
|                                                                                                                    |
| SlepcEigenSolver                                                                                                   |
|   solve_generalized()                  13        12.5278     0.963675    12.5278     0.963675    93.09    93.09    |
 --------------------------------------------------------------------------------------------------------------------
| Totals:                                241394    13.4575                                         100.00            |
 --------------------------------------------------------------------------------------------------------------------

*** Warning, This code is untested, experimental, or likely to see future API changes: src/reduced_basis/rb_parametrized_object.C, line 31, compiled Apr  7 2012 at 15:51:42 ***
mu[0] = 0.2
mu[1] = 0.7
mu[2] = 0.1

output 1, value = 2.35241, bound = 0.00143587
output 2, value = 0.944929, bound = 0.000884036
output 3, value = 0.944929, bound = 0.000884036
output 4, value = 2.35241, bound = 0.00143587

Reading in the basis functions...
Finished reading in the basis functions...

-------------------------------------------------------------------
| Time:           Sat Apr  7 16:01:24 2012                         |
| OS:             Linux                                            |
| HostName:       lkirk-home                                       |
| OS Release:     3.0.0-17-generic                                 |
| OS Version:     #30-Ubuntu SMP Thu Mar 8 20:45:39 UTC 2012       |
| Machine:        x86_64                                           |
| Username:       benkirk                                          |
| Configuration:  ./configure run on Sat Apr  7 15:49:27 CDT 2012  |
-------------------------------------------------------------------
 --------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=0.186489, Active time=0.052959                                               |
 --------------------------------------------------------------------------------------------------------------
| Event                            nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                            w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|--------------------------------------------------------------------------------------------------------------|
|                                                                                                              |
|                                                                                                              |
| DofMap                                                                                                       |
|   add_neighbors_to_send_list()   2         0.0004      0.000183    0.0004      0.000183    0.69     0.69     |
|   compute_sparsity()             2         0.0030      0.001483    0.0039      0.001943    5.60     7.34     |
|   create_dof_constraints()       2         0.0003      0.000159    0.0003      0.000159    0.60     0.60     |
|   distribute_dofs()              2         0.0006      0.000282    0.0022      0.001077    1.06     4.07     |
|   dof_indices()                  3750      0.0020      0.000001    0.0020      0.000001    3.85     3.85     |
|   prepare_send_list()            2         0.0000      0.000001    0.0000      0.000001    0.00     0.00     |
|   reinit()                       2         0.0016      0.000794    0.0016      0.000794    3.00     3.00     |
|                                                                                                              |
| EquationSystems                                                                                              |
|   build_solution_vector()        2         0.0028      0.001412    0.0046      0.002285    5.33     8.63     |
|   read()                         1         0.0186      0.018579    0.0390      0.039025    35.08    73.69    |
|   update()                       1         0.0001      0.000090    0.0001      0.000090    0.17     0.17     |
|                                                                                                              |
| ExodusII_IO                                                                                                  |
|   write_nodal_data()             2         0.0029      0.001466    0.0029      0.001466    5.53     5.53     |
|                                                                                                              |
| Mesh                                                                                                         |
|   find_neighbors()               1         0.0021      0.002093    0.0021      0.002093    3.95     3.95     |
|   renumber_nodes_and_elem()      2         0.0001      0.000066    0.0001      0.000066    0.25     0.25     |
|                                                                                                              |
| MeshCommunication                                                                                            |
|   assign_global_indices()        1         0.0137      0.013670    0.0137      0.013679    25.81    25.83    |
|                                                                                                              |
| MeshOutput                                                                                                   |
|   write_equation_systems()       2         0.0001      0.000026    0.0076      0.003777    0.10     14.27    |
|                                                                                                              |
| MeshTools::Generation                                                                                        |
|   build_cube()                   1         0.0011      0.001120    0.0011      0.001120    2.11     2.11     |
|                                                                                                              |
| Parallel                                                                                                     |
|   allgather()                    26        0.0000      0.000001    0.0000      0.000001    0.05     0.05     |
|   send_receive()                 4         0.0000      0.000001    0.0000      0.000001    0.01     0.01     |
|                                                                                                              |
| Partitioner                                                                                                  |
|   single_partition()             1         0.0001      0.000105    0.0001      0.000105    0.20     0.20     |
|                                                                                                              |
| RBConstruction                                                                                               |
|   clear()                        2         0.0001      0.000058    0.0001      0.000058    0.22     0.22     |
|   load_basis_function()          1         0.0000      0.000024    0.0000      0.000024    0.05     0.05     |
|   load_rb_solution()             1         0.0001      0.000094    0.0001      0.000094    0.18     0.18     |
|                                                                                                              |
| RBEvaluation                                                                                                 |
|   clear()                        1         0.0001      0.000056    0.0001      0.000057    0.11     0.11     |
|   clear_riesz_representors()     2         0.0000      0.000003    0.0000      0.000003    0.01     0.01     |
|   compute_residual_dual_norm()   1         0.0000      0.000033    0.0000      0.000033    0.06     0.06     |
|   rb_solve()                     1         0.0005      0.000506    0.0009      0.000897    0.96     1.69     |
|   read_offline_data_from_files() 1         0.0019      0.001890    0.0020      0.002029    3.57     3.83     |
|   resize_data_structures()       1         0.0001      0.000134    0.0001      0.000139    0.25     0.26     |
|                                                                                                              |
| RBSCMEvaluation                                                                                              |
|   get_SCM_LB()                   1         0.0004      0.000358    0.0004      0.000358    0.68     0.68     |
|   read_offline_data_from_files() 1         0.0003      0.000267    0.0003      0.000267    0.50     0.50     |
 --------------------------------------------------------------------------------------------------------------
| Totals:                          3819      0.0530                                          100.00            |
 --------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running  ./reduced_basis_ex2-opt
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
