<?php $root=""; ?>
<?php require($root."navigation.php"); ?>
<html>
<head>
  <?php load_style($root); ?>
</head>
 
<body>
 
<?php make_navigation("reduced_basis_ex6",$root)?>
 
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
<h1>Reduced Basis: Example 6 - Heat transfer on a curved domain in 3D</h1>


<br><br>In this example we consider heat transfer modeled by a Poisson equation with
Robin boundary condition:
-kappa \Laplacian u = 1, on \Omega
-kappa du\dn = kappa Bi u, on \partial\Omega_Biot,
u = 0 on \partial\Omega_Dirichlet,

<br><br>We consider a reference domain \Omega_hat = [-0.2,0.2]x[-0.2,0.2]x[0,3], and the
physical domain is then obtain via the parametrized mapping:
x = -1/mu + (1/mu+x_hat)*cos(mu*z_hat)
y = y_hat
z = (1/mu+x_hat)*sin(mu*z_hat)
for (x_hat,y_hat,z_hat) \in \Omega_hat. (Here "hats" denotes reference domain.)
Also, the "reference Dirichlet boundaries" are [-0.2,0.2]x[-0.2,0.2]x{0} and
[-0.2,0.2]x[-0.2,0.2]x{3}, and the remaining boundaries are the "Biot" boundaries.


<br><br>Then, after putting the PDE into weak form and mapping it to the reference domain,
we obtain:
\kappa \int_\Omega_hat [ (1+mu*x_hat) v_x w_x + (1+mu*x_hat) v_y w_y + 1/(1+mu*x_hat) v_z w_z ]
+ \kappa Bi \int_\partial\Omega_hat_Biot1 (1-0.2mu) u v
+ \kappa Bi \int_\partial\Omega_hat_Biot2 (1+mu x_hat) u v
+ \kappa Bi \int_\partial\Omega_hat_Biot3 (1+0.2mu) u v
= \int_\Omega_hat (1+mu x_hat) v
where
\partial\Omega_hat_Biot1 = [-0.2] x [-0.2,0.2] x [0,3]
\partial\Omega_hat_Biot2 = [-0.2,0.2] x {-0.2} x [0,3] \UNION [-0.2,0.2] x {0.2} x [0,3]
\partial\Omega_hat_Biot3 = [0.2] x [-0.2,0.2] x [0,3]


<br><br>The term
\kappa \int_\Omega_hat 1/(1+mu*x_hat) v_z w_z 
is "non-affine" (in the Reduced Basis sense), since we can't express it
in the form \sum theta_q(kappa,mu) a(v,w). As a result, (as in
reduced_basis_ex4) we must employ the Empirical Interpolation Method (EIM)
in order to apply the Reduced Basis method here.


<br><br>The approach we use is to construct an EIM approximation, G_EIM, to the vector-valued function
G(x_hat,y_hat;mu) = (1 + mu*x_hat, 1 + mu*x_hat, 1/(1+mu*x_hat))
and then we express the "volumetric integral part" of the left-hand side operator as
a(v,w;mu) = \int_\hat\Omega G_EIM(x_hat,y_hat;mu) \dot (v_x w_x, v_y w_y, v_z w_z).
(We actually only need EIM for the third component of G_EIM, but it's helpful to
demonstrate "vector-valued" EIM here.)


<br><br>

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
        #include "eim_classes.h"
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
Define a function to scale the mesh according to the parameter.
</div>

<div class ="fragment">
<pre>
        void transform_mesh_and_plot(EquationSystems& es, Real curvature, const std::string& filename);
        
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
This is a 3D example
</div>

<div class ="fragment">
<pre>
          libmesh_example_assert(3 == LIBMESH_DIM, "3D support");
          
</pre>
</div>
<div class = "comment">
Parse the input file using GetPot
</div>

<div class ="fragment">
<pre>
          std::string eim_parameters = "eim.in";
          std::string rb_parameters  = "rb.in";
          std::string main_parameters = "reduced_basis_ex6.in";
          GetPot infile(main_parameters);
        
          unsigned int n_elem_xy = infile("n_elem_xy", 1);
          unsigned int n_elem_z  = infile("n_elem_z", 1);
          
</pre>
</div>
<div class = "comment">
Do we write the RB basis functions to disk?
</div>

<div class ="fragment">
<pre>
          bool store_basis_functions = infile("store_basis_functions", true);
        
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
          Mesh mesh;
          MeshTools::Generation::build_cube (mesh,
                                             n_elem_xy, n_elem_xy, n_elem_z,
                                             -0.2, 0.2,
                                             -0.2, 0.2,
                                             0., 3.,
                                             HEX8);
        
</pre>
</div>
<div class = "comment">
Create an equation systems object.
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
          equation_systems.print_info();
          mesh.print_info();
        
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
        
          if(!online_mode) // Perform the Offline stage of the RB method
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
attach the EIM theta objects to the RBEvaluation
</div>

<div class ="fragment">
<pre>
            eim_rb_eval.initialize_eim_theta_objects();
            rb_eval.get_rb_theta_expansion().attach_multiple_A_theta(eim_rb_eval.get_eim_theta_objects());
            
</pre>
</div>
<div class = "comment">
attach the EIM assembly objects to the RBConstruction
</div>

<div class ="fragment">
<pre>
            eim_construction.initialize_eim_assembly_objects();
            rb_construction.get_rb_assembly_expansion().attach_multiple_A_assembly(eim_construction.get_eim_assembly_objects());
        
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
          else // Perform the Online stage of the RB method
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
            rb_eval.get_rb_theta_expansion().attach_multiple_A_theta(eim_rb_eval.get_eim_theta_objects());
            
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
            Real online_curvature = infile("online_curvature", 0.);
            Real online_Bi        = infile("online_Bi", 0.);
            Real online_kappa     = infile("online_kappa", 0.);
            RBParameters online_mu;
            online_mu.set_value("curvature", online_curvature);
            online_mu.set_value("Bi", online_Bi);
            online_mu.set_value("kappa", online_kappa);
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
        
              transform_mesh_and_plot(equation_systems,online_curvature,"RB_sol.e");
            }
          }
        
          return 0;
        }
        
        void transform_mesh_and_plot(EquationSystems& es, Real curvature, const std::string& filename)
        {
</pre>
</div>
<div class = "comment">
Loop over the mesh nodes and move them!
</div>

<div class ="fragment">
<pre>
          MeshBase& mesh = es.get_mesh();
        
          MeshBase::node_iterator       node_it  = mesh.nodes_begin();
          const MeshBase::node_iterator node_end = mesh.nodes_end();
          
          for( ; node_it != node_end; node_it++)
          {
            Node* node = *node_it;
            
            Real x = (*node)(0);
            Real z = (*node)(2);
        
            (*node)(0) = -1./curvature + (1./curvature + x)*cos(curvature*z);
            (*node)(2) = (1./curvature + x)*sin(curvature*z);
          }
        
        #ifdef LIBMESH_HAVE_EXODUS_API
          ExodusII_IO(mesh).write_equation_systems(filename, es);
        #endif
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
  #include <B><FONT COLOR="#BC8F8F">&quot;eim_classes.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;assembly.h&quot;</FONT></B>
  
  using namespace libMesh;
  
  <B><FONT COLOR="#228B22">void</FONT></B> transform_mesh_and_plot(EquationSystems&amp; es, Real curvature, <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp; filename);
  
  <B><FONT COLOR="#228B22">int</FONT></B> main (<B><FONT COLOR="#228B22">int</FONT></B> argc, <B><FONT COLOR="#228B22">char</FONT></B>** argv)
  {
    LibMeshInit init (argc, argv);
  
  #<B><FONT COLOR="#A020F0">if</FONT></B> !defined(LIBMESH_HAVE_XDR)
    libmesh_example_assert(false, <B><FONT COLOR="#BC8F8F">&quot;--enable-xdr&quot;</FONT></B>);
  #elif defined(LIBMESH_DEFAULT_SINGLE_PRECISION)
    libmesh_example_assert(false, <B><FONT COLOR="#BC8F8F">&quot;--disable-singleprecision&quot;</FONT></B>);
  #endif
  
    libmesh_example_assert(3 == LIBMESH_DIM, <B><FONT COLOR="#BC8F8F">&quot;3D support&quot;</FONT></B>);
    
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::string eim_parameters = <B><FONT COLOR="#BC8F8F">&quot;eim.in&quot;</FONT></B>;
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::string rb_parameters  = <B><FONT COLOR="#BC8F8F">&quot;rb.in&quot;</FONT></B>;
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::string main_parameters = <B><FONT COLOR="#BC8F8F">&quot;reduced_basis_ex6.in&quot;</FONT></B>;
    GetPot infile(main_parameters);
  
    <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_elem_xy = infile(<B><FONT COLOR="#BC8F8F">&quot;n_elem_xy&quot;</FONT></B>, 1);
    <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_elem_z  = infile(<B><FONT COLOR="#BC8F8F">&quot;n_elem_z&quot;</FONT></B>, 1);
    
    <B><FONT COLOR="#228B22">bool</FONT></B> store_basis_functions = infile(<B><FONT COLOR="#BC8F8F">&quot;store_basis_functions&quot;</FONT></B>, true);
  
    GetPot command_line (argc, argv);
    <B><FONT COLOR="#228B22">int</FONT></B> online_mode = 0;
    <B><FONT COLOR="#A020F0">if</FONT></B> ( command_line.search(1, <B><FONT COLOR="#BC8F8F">&quot;-online_mode&quot;</FONT></B>) )
      online_mode = command_line.next(online_mode);
  
    Mesh mesh;
    <B><FONT COLOR="#5F9EA0">MeshTools</FONT></B>::Generation::build_cube (mesh,
                                       n_elem_xy, n_elem_xy, n_elem_z,
                                       -0.2, 0.2,
                                       -0.2, 0.2,
                                       0., 3.,
                                       HEX8);
  
    EquationSystems equation_systems (mesh);
  
    SimpleEIMConstruction &amp; eim_construction =
      equation_systems.add_system&lt;SimpleEIMConstruction&gt; (<B><FONT COLOR="#BC8F8F">&quot;EIM&quot;</FONT></B>);
    SimpleRBConstruction &amp; rb_construction =
      equation_systems.add_system&lt;SimpleRBConstruction&gt; (<B><FONT COLOR="#BC8F8F">&quot;RB&quot;</FONT></B>);
  
    equation_systems.init ();
  
    equation_systems.print_info();
    mesh.print_info();
  
    SimpleRBEvaluation rb_eval;
  
    SimpleEIMEvaluation eim_rb_eval;
    
    eim_construction.set_rb_evaluation(eim_rb_eval);
    rb_construction.set_rb_evaluation(rb_eval);
  
    <B><FONT COLOR="#A020F0">if</FONT></B>(!online_mode) <I><FONT COLOR="#B22222">// Perform the Offline stage of the RB method
</FONT></I>    {
      eim_construction.process_parameters_file(eim_parameters);
      eim_construction.print_info();
    
      eim_construction.initialize_rb_construction();
      
      eim_construction.train_reduced_basis();
      eim_construction.get_rb_evaluation().write_offline_data_to_files(<B><FONT COLOR="#BC8F8F">&quot;eim_data&quot;</FONT></B>);
      
      rb_construction.process_parameters_file(rb_parameters);
  
      eim_rb_eval.initialize_eim_theta_objects();
      rb_eval.get_rb_theta_expansion().attach_multiple_A_theta(eim_rb_eval.get_eim_theta_objects());
      
      eim_construction.initialize_eim_assembly_objects();
      rb_construction.get_rb_assembly_expansion().attach_multiple_A_assembly(eim_construction.get_eim_assembly_objects());
  
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
    <B><FONT COLOR="#A020F0">else</FONT></B> <I><FONT COLOR="#B22222">// Perform the Online stage of the RB method
</FONT></I>    {
      eim_rb_eval.read_offline_data_from_files(<B><FONT COLOR="#BC8F8F">&quot;eim_data&quot;</FONT></B>);
  
      eim_rb_eval.initialize_eim_theta_objects();
      rb_eval.get_rb_theta_expansion().attach_multiple_A_theta(eim_rb_eval.get_eim_theta_objects());
      
      rb_eval.read_offline_data_from_files(<B><FONT COLOR="#BC8F8F">&quot;rb_data&quot;</FONT></B>);
  
      Real online_curvature = infile(<B><FONT COLOR="#BC8F8F">&quot;online_curvature&quot;</FONT></B>, 0.);
      Real online_Bi        = infile(<B><FONT COLOR="#BC8F8F">&quot;online_Bi&quot;</FONT></B>, 0.);
      Real online_kappa     = infile(<B><FONT COLOR="#BC8F8F">&quot;online_kappa&quot;</FONT></B>, 0.);
      RBParameters online_mu;
      online_mu.set_value(<B><FONT COLOR="#BC8F8F">&quot;curvature&quot;</FONT></B>, online_curvature);
      online_mu.set_value(<B><FONT COLOR="#BC8F8F">&quot;Bi&quot;</FONT></B>, online_Bi);
      online_mu.set_value(<B><FONT COLOR="#BC8F8F">&quot;kappa&quot;</FONT></B>, online_kappa);
      rb_eval.set_parameters(online_mu);
      rb_eval.print_parameters();
      rb_eval.rb_solve( rb_eval.get_n_basis_functions() );
  
      <B><FONT COLOR="#A020F0">if</FONT></B>(store_basis_functions)
      {
        eim_rb_eval.read_in_basis_functions(eim_construction,<B><FONT COLOR="#BC8F8F">&quot;eim_data&quot;</FONT></B>);
        rb_eval.read_in_basis_functions(rb_construction,<B><FONT COLOR="#BC8F8F">&quot;rb_data&quot;</FONT></B>);
  
        eim_construction.load_rb_solution();
        rb_construction.load_rb_solution();
  
        transform_mesh_and_plot(equation_systems,online_curvature,<B><FONT COLOR="#BC8F8F">&quot;RB_sol.e&quot;</FONT></B>);
      }
    }
  
    <B><FONT COLOR="#A020F0">return</FONT></B> 0;
  }
  
  <B><FONT COLOR="#228B22">void</FONT></B> transform_mesh_and_plot(EquationSystems&amp; es, Real curvature, <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp; filename)
  {
    MeshBase&amp; mesh = es.get_mesh();
  
    <B><FONT COLOR="#5F9EA0">MeshBase</FONT></B>::node_iterator       node_it  = mesh.nodes_begin();
    <B><FONT COLOR="#228B22">const</FONT></B> MeshBase::node_iterator node_end = mesh.nodes_end();
    
    <B><FONT COLOR="#A020F0">for</FONT></B>( ; node_it != node_end; node_it++)
    {
      Node* node = *node_it;
      
      Real x = (*node)(0);
      Real z = (*node)(2);
  
      (*node)(0) = -1./curvature + (1./curvature + x)*cos(curvature*z);
      (*node)(2) = (1./curvature + x)*sin(curvature*z);
    }
  
  #ifdef LIBMESH_HAVE_EXODUS_API
    ExodusII_IO(mesh).write_equation_systems(filename, es);
  #endif
  }
  
</pre> 
<a name="output"></a> 
<br><br><br> <h1> The console output of the program: </h1> 
<pre>
Updated .depend
Compiling C++ (in optimized mode) reduced_basis_ex6.C...
Linking reduced_basis_ex6-opt...
***************************************************************
* Running  ./reduced_basis_ex6-opt
***************************************************************
 
*** Warning, This code is untested, experimental, or likely to see future API changes: src/reduced_basis/rb_parametrized.C, line 40, compiled Aug 22 2012 at 13:57:36 ***
 EquationSystems
  n_systems()=2
   System #0, "EIM"
    Type "RBConstruction"
    Variables="x_comp_of_G" "y_comp_of_G" "z_comp_of_G" 
    Finite Element Types="LAGRANGE" "LAGRANGE" "LAGRANGE" 
    Approximation Orders="FIRST" "FIRST" "FIRST" 
    n_dofs()=18513
    n_local_dofs()=18513
    n_constrained_dofs()=0
    n_local_constrained_dofs()=0
    n_vectors()=1
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 70.545
      Average Off-Processor Bandwidth <= 0
      Maximum  On-Processor Bandwidth <= 81
      Maximum Off-Processor Bandwidth <= 0
    DofMap Constraints
      Number of DoF Constraints = 0
   System #1, "RB"
    Type "RBConstruction"
    Variables="u" 
    Finite Element Types="LAGRANGE" 
    Approximation Orders="FIRST" 
    n_dofs()=6171
    n_local_dofs()=6171
    n_constrained_dofs()=242
    n_local_constrained_dofs()=242
    n_vectors()=1
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 23.515
      Average Off-Processor Bandwidth <= 0
      Maximum  On-Processor Bandwidth <= 27
      Maximum Off-Processor Bandwidth <= 0
    DofMap Constraints
      Number of DoF Constraints = 242
      Average DoF Constraint Length= 0

 Mesh Information:
  mesh_dimension()=3
  spatial_dimension()=3
  n_nodes()=6171
    n_local_nodes()=6171
  n_elem()=5000
    n_local_elem()=5000
    n_active_elem()=5000
  n_subdomains()=1
  n_partitions()=1
  n_processors()=1
  n_threads()=1
  processor_id()=0

Initializing training parameters with deterministic training set...
Parameter curvature: log scaling = 0


RBConstruction parameters:
system name: EIM
constrained_problem: 0
Nmax: 20
Basis training error tolerance: 0.001
Aq operators attached: 0
Fq functions attached: 0
n_outputs: 0
Number of parameters: 1
Parameter curvature: Min = 0.1, Max = 1.0472, value = 0.5
n_training_samples: 25
single-matrix mode? 0
reuse preconditioner? 1
use a relative error bound in greedy? 1
write out data during basis training? 0
quiet mode? 1


RBEIMConstruction parameters:
best fit type: projection

Initializing parametrized functions in training set...
Parametrized functions in training set initialized


---- Performing Greedy basis enrichment ----

---- Basis dimension: 0 ----
Performing truth solve at parameter:
curvature: 1.0472

Enriching the RB space
Updating RB matrices

---- Basis dimension: 1 ----
Performing RB solves on training set
Maximum (relative) error bound is 0.221031

Performing truth solve at parameter:
curvature: 0.1

Enriching the RB space
Updating RB matrices

---- Basis dimension: 2 ----
Performing RB solves on training set
Maximum (relative) error bound is 0.0106172

Performing truth solve at parameter:
curvature: 0.613067

Enriching the RB space
Updating RB matrices

---- Basis dimension: 3 ----
Performing RB solves on training set
Maximum (relative) error bound is 0.000307645

Specified error tolerance reached.
Perform one more Greedy iteration for error bounds.
Performing truth solve at parameter:
curvature: 0.297333

Enriching the RB space
Updating RB matrices

---- Basis dimension: 3 ----
Performing RB solves on training set
Maximum (relative) error bound is 0.000307645

Extra Greedy iteration finished.
Initializing training parameters with random training set...
Parameter Bi: log scaling = 0
Parameter curvature: log scaling = 1
Parameter kappa: log scaling = 1


RBConstruction parameters:
system name: RB
constrained_problem: 0
Nmax: 15
Basis training error tolerance: 0.001
Aq operators attached: 6
Fq functions attached: 2
n_outputs: 0
Number of parameters: 3
Parameter Bi: Min = 0.001, Max = 0.01, value = 0.005
Parameter curvature: Min = 0.1, Max = 1.0472, value = 0.5
Parameter kappa: Min = 0.5, Max = 2, value = 1
n_training_samples: 1000
single-matrix mode? 0
reuse preconditioner? 1
use a relative error bound in greedy? 1
write out data during basis training? 0
quiet mode? 1


---- Performing Greedy basis enrichment ----

---- Basis dimension: 0 ----
Performing RB solves on training set
Maximum (relative) error bound is inf

Performing truth solve at parameter:
Bi: 0.00267569
curvature: 0.53784
kappa: 1.9281

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 1 ----
Performing RB solves on training set
Maximum (relative) error bound is 0.0325359

Performing truth solve at parameter:
Bi: 0.00639032
curvature: 1.04086
kappa: 1.88824

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 2 ----
Performing RB solves on training set
Maximum (relative) error bound is 0.0230728

Performing truth solve at parameter:
Bi: 0.00966468
curvature: 0.177462
kappa: 1.9617

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 3 ----
Performing RB solves on training set
Maximum (relative) error bound is 0.00176716

Performing truth solve at parameter:
Bi: 0.00125625
curvature: 0.133986
kappa: 1.96856

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 4 ----
Performing RB solves on training set
Maximum (relative) error bound is 0.00069616

Specified error tolerance reached.

-------------------------------------------------------------------------------------------------------------------
| Time:           Thu Aug 23 08:06:50 2012                                                                         |
| OS:             Darwin                                                                                           |
| HostName:       www.example.com                                                                                  |
| OS Release:     10.8.0                                                                                           |
| OS Version:     Darwin Kernel Version 10.8.0: Tue Jun  7 16:32:41 PDT 2011; root:xnu-1504.15.3~1/RELEASE_X86_64  |
| Machine:        x86_64                                                                                           |
| Username:       jwpeterson                                                                                       |
| Configuration:  ./configure run on Wed Aug 22 13:54:56 MDT 2012                                                  |
-------------------------------------------------------------------------------------------------------------------
 -------------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=15.5454, Active time=15.2381                                                      |
 -------------------------------------------------------------------------------------------------------------------
| Event                                 nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                                 w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|-------------------------------------------------------------------------------------------------------------------|
|                                                                                                                   |
|                                                                                                                   |
| DofMap                                                                                                            |
|   add_neighbors_to_send_list()        2         0.0030      0.001481    0.0030      0.001481    0.02     0.02     |
|   build_constraint_matrix()           45000     0.0329      0.000001    0.0329      0.000001    0.22     0.22     |
|   build_sparsity()                    2         0.1648      0.082392    0.1736      0.086778    1.08     1.14     |
|   cnstrn_elem_mat_vec()               45000     0.0223      0.000000    0.0223      0.000000    0.15     0.15     |
|   create_dof_constraints()            2         0.0118      0.005897    0.0150      0.007510    0.08     0.10     |
|   distribute_dofs()                   2         0.0052      0.002618    0.0254      0.012702    0.03     0.17     |
|   dof_indices()                       750045    0.5018      0.000001    0.5018      0.000001    3.29     3.29     |
|   prepare_send_list()                 2         0.0000      0.000002    0.0000      0.000002    0.00     0.00     |
|   reinit()                            2         0.0202      0.010082    0.0202      0.010082    0.13     0.13     |
|                                                                                                                   |
| FE                                                                                                                |
|   compute_shape_functions()           434000    2.1375      0.000005    2.1375      0.000005    14.03    14.03    |
|   init_shape_functions()              44078     0.3309      0.000008    0.3309      0.000008    2.17     2.17     |
|   inverse_map()                       360030    0.9691      0.000003    0.9691      0.000003    6.36     6.36     |
|                                                                                                                   |
| FEMap                                                                                                             |
|   compute_affine_map()                434000    0.6465      0.000001    0.6465      0.000001    4.24     4.24     |
|   compute_face_map()                  44000     0.0658      0.000001    0.0658      0.000001    0.43     0.43     |
|   init_face_shape_functions()         20        0.0001      0.000007    0.0001      0.000007    0.00     0.00     |
|   init_reference_to_physical_map()    44078     0.2787      0.000006    0.2787      0.000006    1.83     1.83     |
|                                                                                                                   |
| Mesh                                                                                                              |
|   find_neighbors()                    1         0.0151      0.015077    0.0151      0.015077    0.10     0.10     |
|   renumber_nodes_and_elem()           2         0.0007      0.000334    0.0007      0.000334    0.00     0.00     |
|                                                                                                                   |
| MeshCommunication                                                                                                 |
|   assign_global_indices()             2         0.1583      0.079168    0.1585      0.079227    1.04     1.04     |
|                                                                                                                   |
| MeshTools::Generation                                                                                             |
|   build_cube()                        1         0.0032      0.003185    0.0032      0.003185    0.02     0.02     |
|                                                                                                                   |
| Parallel                                                                                                          |
|   allgather()                         10        0.0000      0.000002    0.0000      0.000002    0.00     0.00     |
|   receive()                           52        0.0001      0.000002    0.0001      0.000002    0.00     0.00     |
|   send()                              52        0.0002      0.000004    0.0002      0.000004    0.00     0.00     |
|   send_receive()                      8         0.0001      0.000013    0.0001      0.000013    0.00     0.00     |
|                                                                                                                   |
| Parallel::Request                                                                                                 |
|   wait()                              52        0.0001      0.000001    0.0001      0.000001    0.00     0.00     |
|                                                                                                                   |
| Partitioner                                                                                                       |
|   single_partition()                  1         0.0003      0.000263    0.0003      0.000263    0.00     0.00     |
|                                                                                                                   |
| PetscLinearSolver                                                                                                 |
|   solve()                             55        2.1190      0.038527    2.1190      0.038527    13.91    13.91    |
|                                                                                                                   |
| PointLocatorTree                                                                                                  |
|   init(no master)                     1         0.0097      0.009747    0.0097      0.009747    0.06     0.06     |
|   operator()                          15        0.0007      0.000048    0.0008      0.000055    0.00     0.01     |
|                                                                                                                   |
| RBConstruction                                                                                                    |
|   add_scaled_matrix_and_vector()      10        1.2950      0.129501    4.1514      0.415143    8.50     27.24    |
|   clear()                             3         0.0094      0.003135    0.0094      0.003135    0.06     0.06     |
|   compute_Fq_representor_innerprods() 2         0.0689      0.034458    0.1667      0.083360    0.45     1.09     |
|   compute_max_error_bound()           9         0.0225      0.002504    1.6863      0.187372    0.15     11.07    |
|   enrich_RB_space()                   4         0.0038      0.000942    0.0038      0.000942    0.02     0.02     |
|   train_reduced_basis()               2         0.0013      0.000660    4.0047      2.002332    0.01     26.28    |
|   truth_assembly()                    4         0.1366      0.034158    0.1367      0.034177    0.90     0.90     |
|   truth_solve()                       4         0.0019      0.000474    0.3195      0.079885    0.01     2.10     |
|   update_RB_system_matrices()         8         0.0640      0.007997    0.0640      0.007997    0.42     0.42     |
|   update_residual_terms()             4         0.1484      0.037104    1.3163      0.329066    0.97     8.64     |
|                                                                                                                   |
| RBEIMConstruction                                                                                                 |
|   compute_best_fit_error()            100       0.6747      0.006747    0.6830      0.006830    4.43     4.48     |
|   enrich_RB_space()                   4         0.1250      0.031262    0.4443      0.111069    0.82     2.92     |
|   truth_solve()                       129       4.1830      0.032427    6.8591      0.053171    27.45    45.01    |
|   update_RB_system_matrices()         4         0.0014      0.000344    0.0302      0.007547    0.01     0.20     |
|                                                                                                                   |
| RBEIMEvaluation                                                                                                   |
|   rb_solve()                          5006      0.0159      0.000003    0.0159      0.000003    0.10     0.10     |
|   write_offline_data_to_files()       1         0.0006      0.000615    0.0194      0.019425    0.00     0.13     |
|                                                                                                                   |
| RBEvaluation                                                                                                      |
|   clear()                             3         0.0001      0.000035    0.0001      0.000035    0.00     0.00     |
|   compute_residual_dual_norm()        5000      0.9190      0.000184    0.9190      0.000184    6.03     6.03     |
|   rb_solve()                          5000      0.0445      0.000009    0.9803      0.000196    0.29     6.43     |
|   resize_data_structures()            2         0.0002      0.000103    0.0002      0.000103    0.00     0.00     |
|   write_offline_data_to_files()       2         0.0236      0.011822    0.0236      0.011822    0.16     0.16     |
 -------------------------------------------------------------------------------------------------------------------
| Totals:                               2215816   15.2381                                         100.00            |
 -------------------------------------------------------------------------------------------------------------------

*** Warning, This code is untested, experimental, or likely to see future API changes: src/reduced_basis/rb_parametrized.C, line 40, compiled Aug 22 2012 at 13:57:36 ***
 EquationSystems
  n_systems()=2
   System #0, "EIM"
    Type "RBConstruction"
    Variables="x_comp_of_G" "y_comp_of_G" "z_comp_of_G" 
    Finite Element Types="LAGRANGE" "LAGRANGE" "LAGRANGE" 
    Approximation Orders="FIRST" "FIRST" "FIRST" 
    n_dofs()=18513
    n_local_dofs()=18513
    n_constrained_dofs()=0
    n_local_constrained_dofs()=0
    n_vectors()=1
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 70.545
      Average Off-Processor Bandwidth <= 0
      Maximum  On-Processor Bandwidth <= 81
      Maximum Off-Processor Bandwidth <= 0
    DofMap Constraints
      Number of DoF Constraints = 0
   System #1, "RB"
    Type "RBConstruction"
    Variables="u" 
    Finite Element Types="LAGRANGE" 
    Approximation Orders="FIRST" 
    n_dofs()=6171
    n_local_dofs()=6171
    n_constrained_dofs()=242
    n_local_constrained_dofs()=242
    n_vectors()=1
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 23.515
      Average Off-Processor Bandwidth <= 0
      Maximum  On-Processor Bandwidth <= 27
      Maximum Off-Processor Bandwidth <= 0
    DofMap Constraints
      Number of DoF Constraints = 242
      Average DoF Constraint Length= 0

 Mesh Information:
  mesh_dimension()=3
  spatial_dimension()=3
  n_nodes()=6171
    n_local_nodes()=6171
  n_elem()=5000
    n_local_elem()=5000
    n_active_elem()=5000
  n_subdomains()=1
  n_partitions()=1
  n_processors()=1
  n_threads()=1
  processor_id()=0

Bi: 0.005
curvature: 1.0472
kappa: 1.3


-------------------------------------------------------------------------------------------------------------------
| Time:           Thu Aug 23 08:06:51 2012                                                                         |
| OS:             Darwin                                                                                           |
| HostName:       inl421321.inl.gov                                                                                |
| OS Release:     10.8.0                                                                                           |
| OS Version:     Darwin Kernel Version 10.8.0: Tue Jun  7 16:32:41 PDT 2011; root:xnu-1504.15.3~1/RELEASE_X86_64  |
| Machine:        x86_64                                                                                           |
| Username:       petejw                                                                                           |
| Configuration:  ./configure run on Wed Aug 22 13:54:56 MDT 2012                                                  |
-------------------------------------------------------------------------------------------------------------------
 --------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=0.484731, Active time=0.434167                                               |
 --------------------------------------------------------------------------------------------------------------
| Event                            nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                            w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|--------------------------------------------------------------------------------------------------------------|
|                                                                                                              |
|                                                                                                              |
| DofMap                                                                                                       |
|   add_neighbors_to_send_list()   2         0.0030      0.001494    0.0030      0.001494    0.69     0.69     |
|   build_sparsity()               2         0.1659      0.082957    0.1748      0.087419    38.21    40.27    |
|   create_dof_constraints()       2         0.0118      0.005888    0.0150      0.007482    2.71     3.45     |
|   distribute_dofs()              2         0.0052      0.002617    0.0254      0.012678    1.21     5.84     |
|   dof_indices()                  35000     0.0207      0.000001    0.0207      0.000001    4.76     4.76     |
|   prepare_send_list()            2         0.0000      0.000002    0.0000      0.000002    0.00     0.00     |
|   reinit()                       2         0.0201      0.010059    0.0201      0.010059    4.63     4.63     |
|                                                                                                              |
| EquationSystems                                                                                              |
|   build_solution_vector()        1         0.0163      0.016289    0.0281      0.028052    3.75     6.46     |
|                                                                                                              |
| ExodusII_IO                                                                                                  |
|   write_nodal_data()             1         0.0078      0.007843    0.0078      0.007843    1.81     1.81     |
|                                                                                                              |
| Mesh                                                                                                         |
|   find_neighbors()               1         0.0151      0.015097    0.0151      0.015097    3.48     3.48     |
|   renumber_nodes_and_elem()      2         0.0007      0.000337    0.0007      0.000337    0.16     0.16     |
|                                                                                                              |
| MeshCommunication                                                                                            |
|   assign_global_indices()        2         0.1599      0.079938    0.1600      0.079997    36.82    36.85    |
|                                                                                                              |
| MeshOutput                                                                                                   |
|   write_equation_systems()       1         0.0001      0.000051    0.0359      0.035946    0.01     8.28     |
|                                                                                                              |
| MeshTools::Generation                                                                                        |
|   build_cube()                   1         0.0032      0.003158    0.0032      0.003158    0.73     0.73     |
|                                                                                                              |
| Parallel                                                                                                     |
|   allgather()                    36        0.0000      0.000001    0.0000      0.000001    0.01     0.01     |
|   send_receive()                 8         0.0001      0.000013    0.0001      0.000013    0.02     0.02     |
|                                                                                                              |
| Partitioner                                                                                                  |
|   single_partition()             1         0.0003      0.000263    0.0003      0.000263    0.06     0.06     |
|                                                                                                              |
| RBConstruction                                                                                               |
|   clear()                        3         0.0024      0.000796    0.0024      0.000796    0.55     0.55     |
|   load_rb_solution()             2         0.0002      0.000094    0.0002      0.000094    0.04     0.04     |
|                                                                                                              |
| RBEIMEvaluation                                                                                              |
|   rb_solve()                     1         0.0001      0.000102    0.0001      0.000102    0.02     0.02     |
|   read_offline_data_from_files() 1         0.0001      0.000115    0.0005      0.000479    0.03     0.11     |
|                                                                                                              |
| RBEvaluation                                                                                                 |
|   clear()                        3         0.0000      0.000013    0.0000      0.000013    0.01     0.01     |
|   compute_residual_dual_norm()   1         0.0005      0.000470    0.0005      0.000470    0.11     0.11     |
|   rb_solve()                     1         0.0000      0.000039    0.0006      0.000611    0.01     0.14     |
|   read_offline_data_from_files() 2         0.0005      0.000271    0.0007      0.000352    0.12     0.16     |
|   resize_data_structures()       2         0.0002      0.000082    0.0002      0.000082    0.04     0.04     |
 --------------------------------------------------------------------------------------------------------------
| Totals:                          35082     0.4342                                          100.00            |
 --------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running  ./reduced_basis_ex6-opt
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
