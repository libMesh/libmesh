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
        
</pre>
</div>
<div class = "comment">
We also need a SCM evaluation object to perform SCM calculations
</div>

<div class ="fragment">
<pre>
          RBSCMEvaluation rb_scm_eval;
          rb_scm_eval.set_rb_theta_expansion( rb_eval.get_rb_theta_expansion() );
        
</pre>
</div>
<div class = "comment">
Tell rb_eval about rb_scm_eval
</div>

<div class ="fragment">
<pre>
          rb_eval.rb_scm_eval = &rb_scm_eval;
          
</pre>
</div>
<div class = "comment">
Finally, need to give rb_scm_con and rb_eval a pointer to the
SCM evaluation object, rb_scm_eval
</div>

<div class ="fragment">
<pre>
          rb_scm_con.set_rb_scm_evaluation(rb_scm_eval);
          
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
            rb_scm_con.process_parameters_file(parameters_filename);
        
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
            rb_con.get_rb_evaluation().write_offline_data_to_files("rb_data");
            rb_scm_con.get_rb_scm_evaluation().write_offline_data_to_files("scm_data");
            
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
              rb_con.get_rb_evaluation().write_out_basis_functions(rb_con,"rb_data");
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
            rb_eval.read_offline_data_from_files("rb_data");
            rb_scm_eval.read_offline_data_from_files("scm_data");
        
</pre>
</div>
<div class = "comment">
Read in online_N and initialize online parameters
</div>

<div class ="fragment">
<pre>
            unsigned int online_N = infile("online_N",1);
            Real online_mu_0 = infile("online_mu_0", 0.);
            Real online_mu_1 = infile("online_mu_1", 0.);
            Real online_mu_2 = infile("online_mu_2", 0.);
            RBParameters online_mu;
            online_mu.set_value("mu_0", online_mu_0);
            online_mu.set_value("mu_1", online_mu_1);
            online_mu.set_value("mu_2", online_mu_2);
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
              rb_eval.read_in_basis_functions(rb_con,"rb_data");
              
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
  
    SimpleRBEvaluation rb_eval;
  
    rb_con.set_rb_evaluation(rb_eval);
  
    RBSCMEvaluation rb_scm_eval;
    rb_scm_eval.set_rb_theta_expansion( rb_eval.get_rb_theta_expansion() );
  
    rb_eval.rb_scm_eval = &amp;rb_scm_eval;
    
    rb_scm_con.set_rb_scm_evaluation(rb_scm_eval);
    
    <B><FONT COLOR="#A020F0">if</FONT></B>(!online_mode) <I><FONT COLOR="#B22222">// Perform the Offline stage of the RB method
</FONT></I>    {
      rb_con.process_parameters_file(parameters_filename);
      rb_scm_con.process_parameters_file(parameters_filename);
  
      rb_con.print_info();
      rb_scm_con.print_info();
  
      rb_con.initialize_rb_construction();
      
      rb_scm_con.perform_SCM_greedy();
  
      rb_con.train_reduced_basis();
      
      rb_con.get_rb_evaluation().write_offline_data_to_files(<B><FONT COLOR="#BC8F8F">&quot;rb_data&quot;</FONT></B>);
      rb_scm_con.get_rb_scm_evaluation().write_offline_data_to_files(<B><FONT COLOR="#BC8F8F">&quot;scm_data&quot;</FONT></B>);
      
      <B><FONT COLOR="#A020F0">if</FONT></B>(store_basis_functions)
      {
        rb_con.get_rb_evaluation().write_out_basis_functions(rb_con,<B><FONT COLOR="#BC8F8F">&quot;rb_data&quot;</FONT></B>);
      }
    }
    <B><FONT COLOR="#A020F0">else</FONT></B> <I><FONT COLOR="#B22222">// Perform the Online stage of the RB method
</FONT></I>    {
  
      rb_eval.read_offline_data_from_files(<B><FONT COLOR="#BC8F8F">&quot;rb_data&quot;</FONT></B>);
      rb_scm_eval.read_offline_data_from_files(<B><FONT COLOR="#BC8F8F">&quot;scm_data&quot;</FONT></B>);
  
      <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> online_N = infile(<B><FONT COLOR="#BC8F8F">&quot;online_N&quot;</FONT></B>,1);
      Real online_mu_0 = infile(<B><FONT COLOR="#BC8F8F">&quot;online_mu_0&quot;</FONT></B>, 0.);
      Real online_mu_1 = infile(<B><FONT COLOR="#BC8F8F">&quot;online_mu_1&quot;</FONT></B>, 0.);
      Real online_mu_2 = infile(<B><FONT COLOR="#BC8F8F">&quot;online_mu_2&quot;</FONT></B>, 0.);
      RBParameters online_mu;
      online_mu.set_value(<B><FONT COLOR="#BC8F8F">&quot;mu_0&quot;</FONT></B>, online_mu_0);
      online_mu.set_value(<B><FONT COLOR="#BC8F8F">&quot;mu_1&quot;</FONT></B>, online_mu_1);
      online_mu.set_value(<B><FONT COLOR="#BC8F8F">&quot;mu_2&quot;</FONT></B>, online_mu_2);
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
        rb_eval.read_in_basis_functions(rb_con,<B><FONT COLOR="#BC8F8F">&quot;rb_data&quot;</FONT></B>);
        
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
Linking reduced_basis_ex2-opt...
***************************************************************
* Running  mpirun -np 6 ./reduced_basis_ex2-opt
***************************************************************
 
*** Warning, This code is untested, experimental, or likely to see future API changes: src/reduced_basis/rb_parametrized.C, line 40, compiled Aug 24 2012 at 15:15:42 ***
 EquationSystems
  n_systems()=2
   System #0, "RBConvectionDiffusion"
    Type "RBConstruction"
    Variables="u" 
    Finite Element Types="LAGRANGE", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" 
    Approximation Orders="FIRST", "THIRD" 
    n_dofs()=676
    n_local_dofs()=129
    n_constrained_dofs()=26
    n_local_constrained_dofs()=11
    n_vectors()=1
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 7.94574
      Average Off-Processor Bandwidth <= 0.550388
      Maximum  On-Processor Bandwidth <= 9
      Maximum Off-Processor Bandwidth <= 5
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
    n_local_dofs()=129
    n_constrained_dofs()=0
    n_local_constrained_dofs()=0
    n_vectors()=0
    n_matrices()=2
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 7.94574
      Average Off-Processor Bandwidth <= 0.550388
      Maximum  On-Processor Bandwidth <= 9
      Maximum Off-Processor Bandwidth <= 5
    DofMap Constraints
      Number of DoF Constraints = 0
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

Initializing training parameters with random training set...
Parameter mu_0: log scaling = 1
Parameter mu_1: log scaling = 1
Parameter mu_2: log scaling = 1

Initializing training parameters with random training set...
Parameter mu_0: log scaling = 1
Parameter mu_1: log scaling = 1
Parameter mu_2: log scaling = 1


RBConstruction parameters:
system name: RBConvectionDiffusion
constrained_problem: 0
Nmax: 10
Basis training error tolerance: 1e-05
Aq operators attached: 3
Fq functions attached: 1
n_outputs: 4
output 0, Q_l = 1
output 1, Q_l = 1
output 2, Q_l = 1
output 3, Q_l = 1
Number of parameters: 3
Parameter mu_0: Min = 0.1, Max = 1, value = 0.2
Parameter mu_1: Min = 0.1, Max = 1, value = 0.7
Parameter mu_2: Min = 0.01, Max = 0.1, value = 0.1
n_training_samples: 100
single-matrix mode? 0
reuse preconditioner? 1
use a relative error bound in greedy? 0
write out data during basis training? 0
quiet mode? 1


RBSCMConstruction parameters:
system name: RBSCMConvectionDiffusion
SCM Greedy tolerance: 0.1
A_q operators attached: 3
Number of parameters: 3
Parameter mu_0: Min = 0.1, Max = 1, value = 0.2
Parameter mu_1: Min = 0.1, Max = 1, value = 0.7
Parameter mu_2: Min = 0.01, Max = 0.1, value = 0.1
n_training_samples: 100


B_min(0) = -1.23204e-15
B_max(0) = 0.999932

B_min(1) = -1.57862e-15
B_max(1) = 0.999933

B_min(2) = -0.380786
B_max(2) = 5.38822e-17

SCM: Added mu = (0.69213,0.521551,0.0434063)

Stability constant for C_J(0) = 0.4222

SCM iteration 0, max_SCM_error = 0.892276

SCM: Added mu = (0.103675,0.923079,0.0210216)

-----------------------------------


Stability constant for C_J(1) = 0.0739609

SCM iteration 1, max_SCM_error = 0.670594

SCM: Added mu = (0.628695,0.103825,0.0939237)

-----------------------------------


Stability constant for C_J(2) = 0.070639

SCM iteration 2, max_SCM_error = 0.345883

SCM: Added mu = (0.100928,0.332112,0.0747578)

-----------------------------------


Stability constant for C_J(3) = 0.0602758

SCM iteration 3, max_SCM_error = 0.180256

SCM: Added mu = (0.188213,0.246457,0.0748834)

-----------------------------------


Stability constant for C_J(4) = 0.118531

SCM iteration 4, max_SCM_error = 0.141211

SCM: Added mu = (0.186783,0.110892,0.083504)

-----------------------------------


Stability constant for C_J(5) = 0.0700704

SCM iteration 5, max_SCM_error = 0.138267

SCM: Added mu = (0.386796,0.775824,0.0109383)

-----------------------------------


Stability constant for C_J(6) = 0.285511

SCM iteration 6, max_SCM_error = 0.115766

SCM: Added mu = (0.182004,0.171219,0.0755275)

-----------------------------------


Stability constant for C_J(7) = 0.10248

SCM iteration 7, max_SCM_error = 0.104045

SCM: Added mu = (0.366062,0.365118,0.0101521)

-----------------------------------


Stability constant for C_J(8) = 0.257414

SCM iteration 8, max_SCM_error = 0.0969933

SCM tolerance of 0.1 reached.

Compute output dual norms
output_dual_innerprods[0][0] = 0.839698
output_dual_innerprods[1][0] = 0.318298
output_dual_innerprods[2][0] = 0.318298
output_dual_innerprods[3][0] = 0.839698

---- Performing Greedy basis enrichment ----

---- Basis dimension: 0 ----
Performing RB solves on training set
Maximum (absolute) error bound is 8.09922

Performing truth solve at parameter:
mu_0: 0.100928
mu_1: 0.332112
mu_2: 0.0747578

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 1 ----
Performing RB solves on training set
Maximum (absolute) error bound is 7.86669

Performing truth solve at parameter:
mu_0: 0.103675
mu_1: 0.923079
mu_2: 0.0210216

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 2 ----
Performing RB solves on training set
Maximum (absolute) error bound is 4.98532

Performing truth solve at parameter:
mu_0: 0.628695
mu_1: 0.103825
mu_2: 0.0939237

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 3 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.646283

Performing truth solve at parameter:
mu_0: 0.825007
mu_1: 0.102381
mu_2: 0.0287278

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 4 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.0354429

Performing truth solve at parameter:
mu_0: 0.166977
mu_1: 0.129849
mu_2: 0.0904392

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 5 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.00783549

Performing truth solve at parameter:
mu_0: 0.878684
mu_1: 0.12689
mu_2: 0.0105435

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 6 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.00581471

Performing truth solve at parameter:
mu_0: 0.15008
mu_1: 0.739223
mu_2: 0.0106897

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 7 ----
Performing RB solves on training set
Maximum (absolute) error bound is 5.89208e-05

Performing truth solve at parameter:
mu_0: 0.136596
mu_1: 0.355732
mu_2: 0.0581777

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 8 ----
Performing RB solves on training set
Maximum (absolute) error bound is 3.11384e-06

Specified error tolerance reached.
In RBSCMEvaluation::write_offline_data_to_files, directory scm_data already exists, overwriting contents.
************************************************************************************************************************
***             WIDEN YOUR WINDOW TO 120 CHARACTERS.  Use 'enscript -r -fCourier9' to print this document            ***
************************************************************************************************************************

---------------------------------------------- PETSc Performance Summary: ----------------------------------------------

./reduced_basis_ex2-opt on a intel-11. named daedalus with 6 processors, by roystgnr Fri Aug 24 15:22:19 2012
Using Petsc Release Version 3.1.0, Patch 5, Mon Sep 27 11:51:54 CDT 2010

                         Max       Max/Min        Avg      Total 
Time (sec):           7.520e+00      1.00022   7.519e+00
Objects:              1.491e+03      1.00000   1.491e+03
Flops:                4.139e+07      1.33176   3.539e+07  2.123e+08
Flops/sec:            5.505e+06      1.33173   4.707e+06  2.824e+07
MPI Messages:         1.825e+04      2.50000   1.096e+04  6.575e+04
MPI Message Lengths:  1.561e+06      2.05077   9.411e+01  6.188e+06
MPI Reductions:       1.274e+04      1.00000

Flop counting convention: 1 flop = 1 real number operation of type (multiply/divide/add/subtract)
                            e.g., VecAXPY() for real vectors of length N --> 2N flops
                            and VecAXPY() for complex vectors of length N --> 8N flops

Summary of Stages:   ----- Time ------  ----- Flops -----  --- Messages ---  -- Message Lengths --  -- Reductions --
                        Avg     %Total     Avg     %Total   counts   %Total     Avg         %Total   counts   %Total 
 0:      Main Stage: 7.5192e+00 100.0%  2.1234e+08 100.0%  6.575e+04 100.0%  9.411e+01      100.0%  1.127e+04  88.4% 

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

STSetUp               15 1.0 1.0362e-03 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
EPSSetUp              15 1.0 2.3120e-01 1.4 0.00e+00 0.0 0.0e+00 0.0e+00 8.0e+02  3  0  0  0  6   3  0  0  0  7     0
EPSSolve              15 1.0 4.8070e+00 1.3 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00 57  0  0  0  0  57  0  0  0  0     0
EPSDense              15 1.0 4.8038e+00 1.3 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00 57  0  0  0  0  57  0  0  0  0     0
VecDot               757 1.0 6.5283e-02 1.8 1.95e+05 1.3 0.0e+00 0.0e+00 7.6e+02  1  0  0  0  6   1  0  0  0  7    16
VecMDot             2619 1.0 1.7087e-01 1.8 9.66e+06 1.3 0.0e+00 0.0e+00 2.6e+03  2 24  0  0 21   2 24  0  0 23   296
VecNorm             2767 1.0 1.3936e-01 1.2 7.14e+05 1.3 0.0e+00 0.0e+00 2.8e+03  2  2  0  0 22   2  2  0  0 25    27
VecScale            2770 1.0 1.3039e-03 1.6 3.56e+05 1.3 0.0e+00 0.0e+00 0.0e+00  0  1  0  0  0   0  1  0  0  0  1432
VecCopy              473 1.0 2.0576e-04 2.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecSet              3911 1.0 1.3561e-03 1.4 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecAXPY              263 1.0 2.3150e-04 1.3 6.79e+04 1.3 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0  1536
VecMAXPY            2730 1.0 4.5249e-03 1.3 1.04e+07 1.3 0.0e+00 0.0e+00 0.0e+00  0 26  0  0  0   0 26  0  0  0 12014
VecAssemblyBegin     250 1.0 7.9370e-02 1.3 0.00e+00 0.0 9.0e+01 1.3e+02 7.5e+02  1  0  0  0  6   1  0  0  0  7     0
VecAssemblyEnd       250 1.0 1.5664e-04 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecScatterBegin     3475 1.0 9.0580e-03 1.8 0.00e+00 0.0 6.3e+04 8.7e+01 0.0e+00  0  0 95 88  0   0  0 95 88  0     0
VecScatterEnd       3475 1.0 1.4493e-01 2.6 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  1  0  0  0  0   1  0  0  0  0     0
VecNormalize        2730 1.0 1.4066e-01 1.2 1.06e+06 1.3 0.0e+00 0.0e+00 2.7e+03  2  3  0  0 21   2  3  0  0 24    39
MatMult             2730 1.0 1.3357e-01 2.5 5.63e+06 1.3 4.9e+04 8.6e+01 0.0e+00  1 14 75 68  0   1 14 75 68  0   222
MatMultAdd           693 1.0 2.4636e-02 1.8 1.52e+06 1.3 1.2e+04 8.6e+01 0.0e+00  0  4 19 17  0   0  4 19 17  0   325
MatSolve            2730 1.0 1.0075e-02 1.4 1.25e+07 1.6 0.0e+00 0.0e+00 0.0e+00  0 29  0  0  0   0 29  0  0  0  6060
MatLUFactorNum        18 1.0 8.6451e-04 1.4 3.63e+05 2.0 0.0e+00 0.0e+00 0.0e+00  0  1  0  0  0   0  1  0  0  0  1904
MatILUFactorSym       18 1.0 1.9913e-03 1.7 0.00e+00 0.0 0.0e+00 0.0e+00 1.8e+01  0  0  0  0  0   0  0  0  0  0     0
MatConvert            30 1.0 4.1169e-02 1.3 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatAssemblyBegin    1094 1.0 2.1644e-01 1.5 0.00e+00 0.0 1.6e+03 4.3e+02 2.0e+03  3  0  3 11 16   3  0  3 11 18     0
MatAssemblyEnd      1094 1.0 8.5563e-02 1.7 0.00e+00 0.0 1.3e+03 2.3e+01 1.2e+03  1  0  2  0 10   1  0  2  0 11     0
MatGetRow          24015 1.0 1.1080e-02 6.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatGetRowIJ           18 1.0 2.5702e-0482.9 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatGetSubMatrice      60 1.0 1.2877e-01 2.3 0.00e+00 0.0 0.0e+00 0.0e+00 2.4e+02  1  0  0  0  2   1  0  0  0  2     0
MatGetOrdering        18 1.0 3.8910e-04 3.0 0.00e+00 0.0 0.0e+00 0.0e+00 3.6e+01  0  0  0  0  0   0  0  0  0  0     0
MatZeroEntries        76 1.0 1.5259e-04 1.3 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPGMRESOrthog      2619 1.0 1.7566e-01 1.7 1.94e+07 1.3 0.0e+00 0.0e+00 2.6e+03  2 48  0  0 21   2 48  0  0 23   577
KSPSetup              70 1.0 2.6178e-04 4.4 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPSolve              37 1.0 3.9519e-01 1.0 3.97e+07 1.3 4.9e+04 8.6e+01 5.4e+03  5 96 75 68 43   5 96 75 68 48   514
PCSetUp               51 1.0 3.9721e-03 1.5 3.63e+05 2.0 0.0e+00 0.0e+00 5.4e+01  0  1  0  0  0   0  1  0  0  0   414
PCSetUpOnBlocks       37 1.0 3.6690e-03 1.5 3.63e+05 2.0 0.0e+00 0.0e+00 5.4e+01  0  1  0  0  0   0  1  0  0  0   449
PCApply             2730 1.0 2.6608e-02 1.3 1.25e+07 1.6 0.0e+00 0.0e+00 0.0e+00  0 29  0  0  0   0 29  0  0  0  2295
------------------------------------------------------------------------------------------------------------------------

Memory usage is given in bytes:

Object Type          Creations   Destructions     Memory  Descendants' Mem.
Reports information only for process 0.

--- Event Stage 0: Main Stage

  Spectral Transform     1              1          568     0
 Eigenproblem Solver     1              1          948     0
       Inner product     1              1          428     0
                 Vec   887            887      1348184     0
         Vec Scatter    41             41        35588     0
           Index Set   286            286       258180     0
   IS L to G Mapping     2              2         2056     0
              Matrix   251            221    204798704     0
       Krylov Solver     3              3        19712     0
      Preconditioner     3              3         2072     0
         PetscRandom    15              1          448     0
========================================================================================================================
Average time to get PetscTime(): 0
Average time for MPI_Barrier(): 4.53949e-05
Average time for zero size MPI_Send(): 8.45194e-05
#PETSc Option Table entries:
-eps_type lapack
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
| Time:           Fri Aug 24 15:22:19 2012                         |
| OS:             Linux                                            |
| HostName:       daedalus                                         |
| OS Release:     2.6.32-34-generic                                |
| OS Version:     #76-Ubuntu SMP Tue Aug 30 17:05:01 UTC 2011      |
| Machine:        x86_64                                           |
| Username:       roystgnr                                         |
| Configuration:  ./configure run on Wed Aug 22 12:44:06 CDT 2012  |
-------------------------------------------------------------------
 --------------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=7.72427, Active time=7.49403                                                       |
 --------------------------------------------------------------------------------------------------------------------
| Event                                  nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                                  w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|--------------------------------------------------------------------------------------------------------------------|
|                                                                                                                    |
|                                                                                                                    |
| CondensedEigenSystem                                                                                               |
|   get_eigenpair()                      15        0.0632      0.004212    2.0488      0.136589    0.84     27.34    |
|   solve()                              15        0.0486      0.003238    4.4924      0.299492    0.65     59.95    |
|                                                                                                                    |
| DofMap                                                                                                             |
|   add_neighbors_to_send_list()         2         0.0002      0.000089    0.0002      0.000116    0.00     0.00     |
|   build_constraint_matrix()            6996      0.0032      0.000000    0.0032      0.000000    0.04     0.04     |
|   build_sparsity()                     2         0.0014      0.000694    0.0019      0.000940    0.02     0.03     |
|   cnstrn_elem_mat_vec()                6996      0.0176      0.000003    0.0176      0.000003    0.23     0.23     |
|   create_dof_constraints()             2         0.0010      0.000488    0.0012      0.000608    0.01     0.02     |
|   distribute_dofs()                    2         0.0005      0.000243    0.0038      0.001880    0.01     0.05     |
|   dof_indices()                        15965     0.0082      0.000001    0.0082      0.000001    0.11     0.11     |
|   prepare_send_list()                  2         0.0000      0.000009    0.0000      0.000009    0.00     0.00     |
|   reinit()                             2         0.0009      0.000461    0.0009      0.000461    0.01     0.01     |
|                                                                                                                    |
| FE                                                                                                                 |
|   compute_shape_functions()            16632     0.0241      0.000001    0.0241      0.000001    0.32     0.32     |
|   init_shape_functions()               2772      0.0089      0.000003    0.0089      0.000003    0.12     0.12     |
|   inverse_map()                        5280      0.0074      0.000001    0.0074      0.000001    0.10     0.10     |
|                                                                                                                    |
| FEMap                                                                                                              |
|   compute_affine_map()                 16632     0.0109      0.000001    0.0109      0.000001    0.15     0.15     |
|   compute_face_map()                   2640      0.0073      0.000003    0.0152      0.000006    0.10     0.20     |
|   init_face_shape_functions()          132       0.0003      0.000003    0.0003      0.000003    0.00     0.00     |
|   init_reference_to_physical_map()     2772      0.0052      0.000002    0.0052      0.000002    0.07     0.07     |
|                                                                                                                    |
| Mesh                                                                                                               |
|   find_neighbors()                     1         0.0005      0.000516    0.0015      0.001480    0.01     0.02     |
|   renumber_nodes_and_elem()            2         0.0000      0.000024    0.0000      0.000024    0.00     0.00     |
|                                                                                                                    |
| MeshCommunication                                                                                                  |
|   assign_global_indices()              1         0.0059      0.005926    0.0120      0.011990    0.08     0.16     |
|   compute_hilbert_indices()            2         0.0023      0.001174    0.0023      0.001174    0.03     0.03     |
|   find_global_indices()                2         0.0003      0.000148    0.0076      0.003779    0.00     0.10     |
|   parallel_sort()                      2         0.0011      0.000548    0.0035      0.001727    0.01     0.05     |
|                                                                                                                    |
| MeshTools::Generation                                                                                              |
|   build_cube()                         1         0.0002      0.000243    0.0002      0.000243    0.00     0.00     |
|                                                                                                                    |
| MetisPartitioner                                                                                                   |
|   partition()                          1         0.0015      0.001492    0.0051      0.005134    0.02     0.07     |
|                                                                                                                    |
| Parallel                                                                                                           |
|   allgather()                          16        0.0034      0.000211    0.0034      0.000211    0.04     0.04     |
|   barrier()                            9         0.0001      0.000010    0.0001      0.000010    0.00     0.00     |
|   broadcast()                          37        0.0007      0.000019    0.0007      0.000019    0.01     0.01     |
|   gather()                             33        0.0039      0.000118    0.0039      0.000118    0.05     0.05     |
|   max(scalar)                          3         0.0018      0.000584    0.0018      0.000584    0.02     0.02     |
|   max(vector)                          3         0.0001      0.000047    0.0001      0.000047    0.00     0.00     |
|   maxloc(scalar)                       18        0.0234      0.001301    0.0234      0.001301    0.31     0.31     |
|   min(vector)                          3         0.0003      0.000088    0.0003      0.000088    0.00     0.00     |
|   probe()                              90        0.0063      0.000069    0.0063      0.000069    0.08     0.08     |
|   receive()                            282       0.0002      0.000001    0.0065      0.000023    0.00     0.09     |
|   send()                               122       0.0002      0.000001    0.0002      0.000001    0.00     0.00     |
|   send_receive()                       98        0.0002      0.000002    0.0068      0.000069    0.00     0.09     |
|   sum()                                47        1.9924      0.042392    1.9924      0.042392    26.59    26.59    |
|                                                                                                                    |
| Parallel::Request                                                                                                  |
|   wait()                               282       0.0002      0.000001    0.0002      0.000001    0.00     0.00     |
|                                                                                                                    |
| Partitioner                                                                                                        |
|   set_node_processor_ids()             1         0.0001      0.000095    0.0011      0.001079    0.00     0.01     |
|   set_parent_processor_ids()           1         0.0000      0.000040    0.0000      0.000040    0.00     0.00     |
|                                                                                                                    |
| PetscLinearSolver                                                                                                  |
|   solve()                              37        0.4110      0.011107    0.4110      0.011107    5.48     5.48     |
|                                                                                                                    |
| RBConstruction                                                                                                     |
|   add_scaled_Aq()                      57        0.0010      0.000017    0.2494      0.004376    0.01     3.33     |
|   add_scaled_matrix_and_vector()       66        0.2006      0.003039    0.2989      0.004529    2.68     3.99     |
|   clear()                              1         0.0003      0.000290    0.0003      0.000290    0.00     0.00     |
|   compute_Fq_representor_innerprods()  1         0.0007      0.000729    0.0107      0.010723    0.01     0.14     |
|   compute_max_error_bound()            9         0.0008      0.000091    0.0339      0.003765    0.01     0.45     |
|   compute_output_dual_innerprods()     1         0.0037      0.003689    0.0500      0.050028    0.05     0.67     |
|   enrich_RB_space()                    8         0.0110      0.001370    0.0110      0.001370    0.15     0.15     |
|   train_reduced_basis()                1         0.0016      0.001562    0.5867      0.586673    0.02     7.83     |
|   truth_assembly()                     8         0.0098      0.001231    0.0098      0.001231    0.13     0.13     |
|   truth_solve()                        8         0.0046      0.000580    0.1183      0.014784    0.06     1.58     |
|   update_RB_system_matrices()          8         0.0335      0.004183    0.0335      0.004183    0.45     0.45     |
|   update_residual_terms()              8         0.0769      0.009612    0.3277      0.040967    1.03     4.37     |
|                                                                                                                    |
| RBEvaluation                                                                                                       |
|   clear()                              1         0.0000      0.000031    0.0000      0.000031    0.00     0.00     |
|   compute_residual_dual_norm()         153       0.0057      0.000038    0.0057      0.000038    0.08     0.08     |
|   rb_solve()                           153       0.0017      0.000011    0.0137      0.000090    0.02     0.18     |
|   resize_data_structures()             1         0.0000      0.000038    0.0000      0.000038    0.00     0.00     |
|   write_offline_data_to_files()        1         0.0011      0.001065    0.0011      0.001065    0.01     0.01     |
|                                                                                                                    |
| RBSCMConstruction                                                                                                  |
|   add_scaled_symm_Aq()                 57        0.0001      0.000002    0.2495      0.004378    0.00     3.33     |
|   compute_SCM_bounding_box()           1         0.0007      0.000694    2.3715      2.371532    0.01     31.65    |
|   compute_SCM_bounds_on_training_set() 9         0.0009      0.000098    0.0141      0.001572    0.01     0.19     |
|   enrich_C_J()                         9         0.0040      0.000443    0.0044      0.000494    0.05     0.06     |
|   evaluate_stability_constant()        9         0.0102      0.001130    4.4301      0.492238    0.14     59.12    |
|   perform_SCM_greedy()                 1         0.0044      0.004439    6.8247      6.824714    0.06     91.07    |
|                                                                                                                    |
| RBSCMEvaluation                                                                                                    |
|   get_SCM_LB()                         306       0.0112      0.000037    0.0112      0.000037    0.15     0.15     |
|   get_SCM_UB()                         153       0.0004      0.000003    0.0004      0.000003    0.01     0.01     |
|   write_offline_data_to_files()        1         0.0003      0.000277    0.0003      0.000277    0.00     0.00     |
|                                                                                                                    |
| SlepcEigenSolver                                                                                                   |
|   solve_generalized()                  15        4.4438      0.296254    4.4438      0.296254    59.30    59.30    |
 --------------------------------------------------------------------------------------------------------------------
| Totals:                                78999     7.4940                                          100.00            |
 --------------------------------------------------------------------------------------------------------------------

*** Warning, This code is untested, experimental, or likely to see future API changes: src/reduced_basis/rb_parametrized.C, line 40, compiled Aug 24 2012 at 15:15:42 ***
 EquationSystems
  n_systems()=2
   System #0, "RBConvectionDiffusion"
    Type "RBConstruction"
    Variables="u" 
    Finite Element Types="LAGRANGE", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" 
    Approximation Orders="FIRST", "THIRD" 
    n_dofs()=676
    n_local_dofs()=129
    n_constrained_dofs()=26
    n_local_constrained_dofs()=11
    n_vectors()=1
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 7.94574
      Average Off-Processor Bandwidth <= 0.550388
      Maximum  On-Processor Bandwidth <= 9
      Maximum Off-Processor Bandwidth <= 5
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
    n_local_dofs()=129
    n_constrained_dofs()=0
    n_local_constrained_dofs()=0
    n_vectors()=0
    n_matrices()=2
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 7.94574
      Average Off-Processor Bandwidth <= 0.550388
      Maximum  On-Processor Bandwidth <= 9
      Maximum Off-Processor Bandwidth <= 5
    DofMap Constraints
      Number of DoF Constraints = 0
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

mu_0: 0.2
mu_1: 0.7
mu_2: 0.1

output 1, value = 2.35244, bound = 0.00333093
output 2, value = 0.944998, bound = 0.00205079
output 3, value = 0.944998, bound = 0.00205079
output 4, value = 2.35244, bound = 0.00333093

************************************************************************************************************************
***             WIDEN YOUR WINDOW TO 120 CHARACTERS.  Use 'enscript -r -fCourier9' to print this document            ***
************************************************************************************************************************

---------------------------------------------- PETSc Performance Summary: ----------------------------------------------

./reduced_basis_ex2-opt on a intel-11. named daedalus with 6 processors, by roystgnr Fri Aug 24 15:22:19 2012
Using Petsc Release Version 3.1.0, Patch 5, Mon Sep 27 11:51:54 CDT 2010

                         Max       Max/Min        Avg      Total 
Time (sec):           1.055e-01      1.00403   1.051e-01
Objects:              3.400e+01      1.00000   3.400e+01
Flops:                1.548e+03      1.32990   1.352e+03  8.112e+03
Flops/sec:            1.468e+04      1.32456   1.286e+04  7.715e+04
MPI Messages:         3.000e+01      2.50000   1.900e+01  1.140e+02
MPI Message Lengths:  1.840e+03      1.95745   6.449e+01  7.352e+03
MPI Reductions:       8.400e+01      1.00000

Flop counting convention: 1 flop = 1 real number operation of type (multiply/divide/add/subtract)
                            e.g., VecAXPY() for real vectors of length N --> 2N flops
                            and VecAXPY() for complex vectors of length N --> 8N flops

Summary of Stages:   ----- Time ------  ----- Flops -----  --- Messages ---  -- Message Lengths --  -- Reductions --
                        Avg     %Total     Avg     %Total   counts   %Total     Avg         %Total   counts   %Total 
 0:      Main Stage: 1.0508e-01  99.9%  8.1120e+03 100.0%  1.140e+02 100.0%  6.449e+01      100.0%  4.900e+01  58.3% 

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

VecCopy                3 1.0 5.9605e-06 6.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecSet                16 1.0 1.4067e-05 2.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecAXPY                6 1.0 5.3167e-05 2.0 1.55e+03 1.3 0.0e+00 0.0e+00 0.0e+00  0100  0  0  0   0100  0  0  0   153
VecAssemblyBegin      11 1.0 2.3223e-0217.6 0.00e+00 0.0 0.0e+00 0.0e+00 3.3e+01 18  0  0  0 39  18  0  0  0 67     0
VecAssemblyEnd        11 1.0 2.9087e-05 1.4 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecScatterBegin        2 1.0 3.6955e-05 1.4 0.00e+00 0.0 3.8e+01 1.3e+02 0.0e+00  0  0 33 65  0   0  0 33 65  0     0
VecScatterEnd          2 1.0 2.1935e-05 3.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatZeroEntries         6 1.0 1.7881e-05 1.5 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
------------------------------------------------------------------------------------------------------------------------

Memory usage is given in bytes:

Object Type          Creations   Destructions     Memory  Descendants' Mem.
Reports information only for process 0.

--- Event Stage 0: Main Stage

                 Vec    17             17        38632     0
         Vec Scatter     2              2         1736     0
           Index Set     4              4         2280     0
   IS L to G Mapping     2              2         2056     0
              Matrix     9              9        68352     0
========================================================================================================================
Average time to get PetscTime(): 0
Average time for MPI_Barrier(): 4.35829e-05
Average time for zero size MPI_Send(): 4.28359e-05
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
| Time:           Fri Aug 24 15:22:19 2012                         |
| OS:             Linux                                            |
| HostName:       daedalus                                         |
| OS Release:     2.6.32-34-generic                                |
| OS Version:     #76-Ubuntu SMP Tue Aug 30 17:05:01 UTC 2011      |
| Machine:        x86_64                                           |
| Username:       roystgnr                                         |
| Configuration:  ./configure run on Wed Aug 22 12:44:06 CDT 2012  |
-------------------------------------------------------------------
 --------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=0.318532, Active time=0.071575                                               |
 --------------------------------------------------------------------------------------------------------------
| Event                            nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                            w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|--------------------------------------------------------------------------------------------------------------|
|                                                                                                              |
|                                                                                                              |
| DofMap                                                                                                       |
|   add_neighbors_to_send_list()   2         0.0002      0.000076    0.0002      0.000095    0.21     0.27     |
|   build_sparsity()               2         0.0023      0.001153    0.0029      0.001430    3.22     4.00     |
|   create_dof_constraints()       2         0.0021      0.001027    0.0023      0.001131    2.87     3.16     |
|   distribute_dofs()              2         0.0004      0.000200    0.0029      0.001471    0.56     4.11     |
|   dof_indices()                  2397      0.0007      0.000000    0.0007      0.000000    0.98     0.98     |
|   prepare_send_list()            2         0.0000      0.000003    0.0000      0.000003    0.01     0.01     |
|   reinit()                       2         0.0008      0.000398    0.0008      0.000398    1.11     1.11     |
|                                                                                                              |
| EquationSystems                                                                                              |
|   build_solution_vector()        2         0.0006      0.000316    0.0018      0.000906    0.88     2.53     |
|                                                                                                              |
| ExodusII_IO                                                                                                  |
|   write_nodal_data()             2         0.0019      0.000954    0.0019      0.000954    2.67     2.67     |
|                                                                                                              |
| Mesh                                                                                                         |
|   find_neighbors()               1         0.0005      0.000508    0.0017      0.001685    0.71     2.35     |
|   renumber_nodes_and_elem()      2         0.0001      0.000026    0.0001      0.000026    0.07     0.07     |
|                                                                                                              |
| MeshCommunication                                                                                            |
|   assign_global_indices()        1         0.0317      0.031716    0.0346      0.034559    44.31    48.28    |
|   compute_hilbert_indices()      2         0.0023      0.001138    0.0023      0.001138    3.18     3.18     |
|   find_global_indices()          2         0.0002      0.000121    0.0072      0.003618    0.34     10.11    |
|   parallel_sort()                2         0.0013      0.000649    0.0039      0.001932    1.81     5.40     |
|                                                                                                              |
| MeshOutput                                                                                                   |
|   write_equation_systems()       2         0.0000      0.000011    0.0037      0.001871    0.03     5.23     |
|                                                                                                              |
| MeshTools::Generation                                                                                        |
|   build_cube()                   1         0.0002      0.000243    0.0002      0.000243    0.34     0.34     |
|                                                                                                              |
| MetisPartitioner                                                                                             |
|   partition()                    1         0.0028      0.002831    0.0063      0.006300    3.96     8.80     |
|                                                                                                              |
| Parallel                                                                                                     |
|   allgather()                    32        0.0025      0.000079    0.0025      0.000079    3.52     3.52     |
|   barrier()                      1         0.0013      0.001311    0.0013      0.001311    1.83     1.83     |
|   broadcast()                    134       0.0007      0.000005    0.0006      0.000005    0.93     0.90     |
|   gather()                       1         0.0000      0.000005    0.0000      0.000005    0.01     0.01     |
|   max(scalar)                    3         0.0014      0.000482    0.0014      0.000482    2.02     2.02     |
|   max(vector)                    3         0.0001      0.000041    0.0001      0.000041    0.17     0.17     |
|   min(vector)                    3         0.0020      0.000683    0.0020      0.000683    2.86     2.86     |
|   probe()                        90        0.0011      0.000012    0.0011      0.000012    1.50     1.50     |
|   receive()                      90        0.0008      0.000009    0.0019      0.000021    1.17     2.68     |
|   send()                         90        0.0001      0.000001    0.0001      0.000001    0.10     0.10     |
|   send_receive()                 98        0.0004      0.000004    0.0025      0.000025    0.58     3.44     |
|   sum()                          22        0.0038      0.000172    0.0038      0.000172    5.29     5.29     |
|                                                                                                              |
| Parallel::Request                                                                                            |
|   wait()                         90        0.0000      0.000000    0.0000      0.000000    0.06     0.06     |
|                                                                                                              |
| Partitioner                                                                                                  |
|   set_node_processor_ids()       1         0.0001      0.000082    0.0002      0.000186    0.11     0.26     |
|   set_parent_processor_ids()     1         0.0000      0.000040    0.0000      0.000040    0.06     0.06     |
|                                                                                                              |
| RBConstruction                                                                                               |
|   clear()                        1         0.0001      0.000097    0.0001      0.000097    0.14     0.14     |
|   load_basis_function()          1         0.0000      0.000045    0.0000      0.000045    0.06     0.06     |
|   load_rb_solution()             1         0.0003      0.000341    0.0003      0.000341    0.48     0.48     |
|                                                                                                              |
| RBEvaluation                                                                                                 |
|   clear()                        1         0.0000      0.000043    0.0000      0.000043    0.06     0.06     |
|   compute_residual_dual_norm()   1         0.0001      0.000061    0.0001      0.000061    0.09     0.09     |
|   rb_solve()                     1         0.0069      0.006868    0.0071      0.007128    9.60     9.96     |
|   read_offline_data_from_files() 1         0.0012      0.001196    0.0012      0.001235    1.67     1.73     |
|   resize_data_structures()       1         0.0000      0.000039    0.0000      0.000039    0.05     0.05     |
|                                                                                                              |
| RBSCMEvaluation                                                                                              |
|   get_SCM_LB()                   1         0.0002      0.000199    0.0002      0.000199    0.28     0.28     |
|   read_offline_data_from_files() 1         0.0001      0.000081    0.0001      0.000081    0.11     0.11     |
 --------------------------------------------------------------------------------------------------------------
| Totals:                          3096      0.0716                                          100.00            |
 --------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running  mpirun -np 6 ./reduced_basis_ex2-opt
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
