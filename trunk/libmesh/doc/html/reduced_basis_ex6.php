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
Linking reduced_basis_ex6-opt...
***************************************************************
* Running  mpirun -np 6 ./reduced_basis_ex6-opt -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
***************************************************************
 
*** Warning, This code is untested, experimental, or likely to see future API changes: src/reduced_basis/rb_parametrized.C, line 40, compiled Aug 24 2012 at 15:15:42 ***
 EquationSystems
  n_systems()=2
   System #0, "EIM"
    Type "RBConstruction"
    Variables="x_comp_of_G" "y_comp_of_G" "z_comp_of_G" 
    Finite Element Types="LAGRANGE", "JACOBI_20_00" "LAGRANGE", "JACOBI_20_00" "LAGRANGE", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" "CARTESIAN" "CARTESIAN" 
    Approximation Orders="FIRST", "THIRD" "FIRST", "THIRD" "FIRST", "THIRD" 
    n_dofs()=18513
    n_local_dofs()=3414
    n_constrained_dofs()=0
    n_local_constrained_dofs()=0
    n_vectors()=1
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 22.0896
      Average Off-Processor Bandwidth <= 1.72056
      Maximum  On-Processor Bandwidth <= 27
      Maximum Off-Processor Bandwidth <= 14
    DofMap Constraints
      Number of DoF Constraints = 0
      Number of Node Constraints = 0
   System #1, "RB"
    Type "RBConstruction"
    Variables="u" 
    Finite Element Types="LAGRANGE", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" 
    Approximation Orders="FIRST", "THIRD" 
    n_dofs()=6171
    n_local_dofs()=1138
    n_constrained_dofs()=243
    n_local_constrained_dofs()=1
    n_vectors()=1
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 22.0896
      Average Off-Processor Bandwidth <= 1.72056
      Maximum  On-Processor Bandwidth <= 27
      Maximum Off-Processor Bandwidth <= 14
    DofMap Constraints
      Number of DoF Constraints = 242
      Average DoF Constraint Length= 0
      Number of Node Constraints = 0

 Mesh Information:
  mesh_dimension()=3
  spatial_dimension()=3
  n_nodes()=6171
    n_local_nodes()=1138
  n_elem()=5000
    n_local_elem()=833
    n_active_elem()=5000
  n_subdomains()=1
  n_partitions()=6
  n_processors()=6
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
Bi: 0.00504949
curvature: 0.312457
kappa: 1.2467

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 1 ----
Performing RB solves on training set
Maximum (relative) error bound is 0.0403079

Performing truth solve at parameter:
Bi: 0.00947927
curvature: 0.919606
kappa: 1.97243

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 2 ----
Performing RB solves on training set
Maximum (relative) error bound is 0.0162604

Performing truth solve at parameter:
Bi: 0.00166235
curvature: 0.854531
kappa: 1.82214

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 3 ----
Performing RB solves on training set
Maximum (relative) error bound is 0.000971329

Specified error tolerance reached.
************************************************************************************************************************
***             WIDEN YOUR WINDOW TO 120 CHARACTERS.  Use 'enscript -r -fCourier9' to print this document            ***
************************************************************************************************************************

---------------------------------------------- PETSc Performance Summary: ----------------------------------------------

./reduced_basis_ex6-opt on a intel-11. named daedalus with 6 processors, by roystgnr Fri Aug 24 15:26:09 2012
Using Petsc Release Version 3.1.0, Patch 5, Mon Sep 27 11:51:54 CDT 2010

                         Max       Max/Min        Avg      Total 
Time (sec):           9.724e+00      1.00024   9.722e+00
Objects:              3.280e+02      1.00000   3.280e+02
Flops:                1.364e+09      1.39770   1.208e+09  7.246e+09
Flops/sec:            1.403e+08      1.39769   1.242e+08  7.453e+08
MPI Messages:         5.809e+03      1.93504   4.873e+03  2.924e+04
MPI Message Lengths:  1.366e+07      1.76179   2.319e+03  6.781e+07
MPI Reductions:       8.308e+03      1.00000

Flop counting convention: 1 flop = 1 real number operation of type (multiply/divide/add/subtract)
                            e.g., VecAXPY() for real vectors of length N --> 2N flops
                            and VecAXPY() for complex vectors of length N --> 8N flops

Summary of Stages:   ----- Time ------  ----- Flops -----  --- Messages ---  -- Message Lengths --  -- Reductions --
                        Avg     %Total     Avg     %Total   counts   %Total     Avg         %Total   counts   %Total 
 0:      Main Stage: 9.7219e+00 100.0%  7.2457e+09 100.0%  2.924e+04 100.0%  2.319e+03      100.0%  8.138e+03  98.0% 

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

VecDot               531 1.0 7.3075e-02 1.9 2.30e+06 1.3 0.0e+00 0.0e+00 5.3e+02  1  0  0  0  6   1  0  0  0  7   169
VecMDot             1973 1.0 9.2566e-01 1.6 7.15e+07 1.3 0.0e+00 0.0e+00 2.0e+03  8  5  0  0 24   8  5  0  0 24   415
VecNorm             2215 1.0 1.0698e+00 1.8 6.73e+06 1.3 0.0e+00 0.0e+00 2.2e+03  8  0  0  0 27   8  0  0  0 27    34
VecScale            2098 1.0 3.4974e-03 1.6 3.29e+06 1.3 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0  5052
VecCopy              605 1.0 2.0566e-03 1.3 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecSet              2866 1.0 5.1210e-03 1.3 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecAXPY              430 1.0 1.5103e-02 2.0 2.28e+06 1.3 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0   811
VecMAXPY            2067 1.0 3.4005e-02 1.4 7.77e+07 1.3 0.0e+00 0.0e+00 0.0e+00  0  6  0  0  0   0  6  0  0  0 12263
VecAssemblyBegin     427 1.0 3.2549e-01 2.6 0.00e+00 0.0 2.7e+02 7.8e+03 1.3e+03  2  0  1  3 15   2  0  1  3 16     0
VecAssemblyEnd       427 1.0 5.8389e-04 1.7 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecScatterBegin     2751 1.0 1.9445e-02 1.8 0.00e+00 0.0 2.8e+04 2.1e+03 0.0e+00  0  0 95 85  0   0  0 95 85  0     0
VecScatterEnd       2751 1.0 6.8691e-0117.5 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  4  0  0  0  0   4  0  0  0  0     0
VecNormalize        2067 1.0 9.7514e-01 1.8 9.75e+06 1.3 0.0e+00 0.0e+00 2.1e+03  7  1  0  0 25   7  1  0  0 25    54
MatMult             2067 1.0 8.5584e-01 4.3 1.50e+08 1.3 2.1e+04 1.5e+03 0.0e+00  6 11 71 45  0   6 11 71 45  0   938
MatMultAdd           512 1.0 6.3750e-02 1.8 5.31e+07 1.3 5.1e+03 2.1e+03 0.0e+00  1  4 18 16  0   1  4 18 16  0  4461
MatSolve            2067 1.0 1.0982e+00 1.4 8.80e+08 1.4 0.0e+00 0.0e+00 0.0e+00  9 64  0  0  0   9 64  0  0  0  4254
MatLUFactorNum         8 1.0 1.3310e-01 2.1 1.19e+08 1.5 0.0e+00 0.0e+00 0.0e+00  1  8  0  0  0   1  8  0  0  0  4561
MatILUFactorSym        8 1.0 3.7151e-01 1.8 0.00e+00 0.0 0.0e+00 0.0e+00 8.0e+00  3  0  0  0  0   3  0  0  0  0     0
MatAssemblyBegin     651 1.0 2.1279e+00 2.3 0.00e+00 0.0 1.2e+02 2.4e+04 1.3e+03 16  0  0  4 16  16  0  0  4 16     0
MatAssemblyEnd       651 1.0 1.9915e-01 1.6 0.00e+00 0.0 2.0e+02 3.8e+02 7.2e+02  2  0  1  0  9   2  0  1  0  9     0
MatGetRow          32200 1.3 1.6208e-02 3.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatGetRowIJ            8 1.0 1.2875e-05 4.5 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatGetOrdering         8 1.0 1.6999e-04 1.3 0.00e+00 0.0 0.0e+00 0.0e+00 1.6e+01  0  0  0  0  0   0  0  0  0  0     0
MatZeroEntries        37 1.0 2.2850e-03 1.4 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPGMRESOrthog      1973 1.0 9.5648e-01 1.6 1.43e+08 1.3 0.0e+00 0.0e+00 2.0e+03  8 11  0  0 24   8 11  0  0 24   803
KSPSetup              56 1.0 3.6955e-04 4.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPSolve              48 1.0 3.4485e+00 1.0 1.31e+09 1.4 2.1e+04 1.5e+03 4.1e+03 35 96 71 45 49  35 96 71 45 51  2013
PCSetUp               16 1.0 4.9368e-01 1.8 1.19e+08 1.5 0.0e+00 0.0e+00 2.4e+01  4  8  0  0  0   4  8  0  0  0  1230
PCSetUpOnBlocks       48 1.0 4.9331e-01 1.8 1.19e+08 1.5 0.0e+00 0.0e+00 2.4e+01  4  8  0  0  0   4  8  0  0  0  1231
PCApply             2067 1.0 1.1241e+00 1.4 8.80e+08 1.4 0.0e+00 0.0e+00 0.0e+00  9 64  0  0  0   9 64  0  0  0  4156
------------------------------------------------------------------------------------------------------------------------

Memory usage is given in bytes:

Object Type          Creations   Destructions     Memory  Descendants' Mem.
Reports information only for process 0.

--- Event Stage 0: Main Stage

                 Vec   174            174      2843448     0
         Vec Scatter    30             30        26040     0
           Index Set    71             71       201088     0
   IS L to G Mapping     5              5        74508     0
              Matrix    40             40     23972056     0
       Krylov Solver     4              4        37760     0
      Preconditioner     4              4         2816     0
========================================================================================================================
Average time to get PetscTime(): 0
Average time for MPI_Barrier(): 4.18186e-05
Average time for zero size MPI_Send(): 4.84784e-05
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
| Time:           Fri Aug 24 15:26:09 2012                         |
| OS:             Linux                                            |
| HostName:       daedalus                                         |
| OS Release:     2.6.32-34-generic                                |
| OS Version:     #76-Ubuntu SMP Tue Aug 30 17:05:01 UTC 2011      |
| Machine:        x86_64                                           |
| Username:       roystgnr                                         |
| Configuration:  ./configure run on Wed Aug 22 12:44:06 CDT 2012  |
-------------------------------------------------------------------
 -------------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=9.93006, Active time=7.92255                                                      |
 -------------------------------------------------------------------------------------------------------------------
| Event                                 nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                                 w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|-------------------------------------------------------------------------------------------------------------------|
|                                                                                                                   |
|                                                                                                                   |
| DofMap                                                                                                            |
|   add_neighbors_to_send_list()        2         0.0016      0.000796    0.0021      0.001046    0.02     0.03     |
|   build_constraint_matrix()           7497      0.0026      0.000000    0.0026      0.000000    0.03     0.03     |
|   build_sparsity()                    2         0.0277      0.013826    0.0360      0.017988    0.35     0.45     |
|   cnstrn_elem_mat_vec()               7497      0.0013      0.000000    0.0013      0.000000    0.02     0.02     |
|   create_dof_constraints()            2         0.0085      0.004239    0.0104      0.005212    0.11     0.13     |
|   distribute_dofs()                   2         0.0047      0.002361    0.0688      0.034419    0.06     0.87     |
|   dof_indices()                       148348    0.0681      0.000000    0.0681      0.000000    0.86     0.86     |
|   prepare_send_list()                 2         0.0002      0.000120    0.0002      0.000120    0.00     0.00     |
|   reinit()                            2         0.0129      0.006441    0.0129      0.006441    0.16     0.16     |
|                                                                                                                   |
| FE                                                                                                                |
|   compute_shape_functions()           71714     0.4438      0.000006    0.4438      0.000006    5.60     5.60     |
|   init_shape_functions()              6818      0.0669      0.000010    0.0669      0.000010    0.84     0.84     |
|   inverse_map()                       60006     0.1328      0.000002    0.1328      0.000002    1.68     1.68     |
|                                                                                                                   |
| FEMap                                                                                                             |
|   compute_affine_map()                71714     0.0983      0.000001    0.0983      0.000001    1.24     1.24     |
|   compute_face_map()                  6740      0.0093      0.000001    0.0093      0.000001    0.12     0.12     |
|   init_face_shape_functions()         20        0.0001      0.000007    0.0001      0.000007    0.00     0.00     |
|   init_reference_to_physical_map()    6818      0.0434      0.000006    0.0434      0.000006    0.55     0.55     |
|                                                                                                                   |
| Mesh                                                                                                              |
|   find_neighbors()                    1         0.0081      0.008096    0.0190      0.019015    0.10     0.24     |
|   renumber_nodes_and_elem()           2         0.0007      0.000371    0.0007      0.000371    0.01     0.01     |
|                                                                                                                   |
| MeshCommunication                                                                                                 |
|   assign_global_indices()             2         0.0883      0.044168    0.1688      0.084386    1.11     2.13     |
|   compute_hilbert_indices()           2         0.0179      0.008961    0.0179      0.008961    0.23     0.23     |
|   find_global_indices()               2         0.0019      0.000942    0.0417      0.020847    0.02     0.53     |
|   parallel_sort()                     2         0.0017      0.000825    0.0198      0.009921    0.02     0.25     |
|                                                                                                                   |
| MeshTools::Generation                                                                                             |
|   build_cube()                        1         0.0022      0.002216    0.0022      0.002216    0.03     0.03     |
|                                                                                                                   |
| MetisPartitioner                                                                                                  |
|   partition()                         1         0.0203      0.020291    0.0414      0.041391    0.26     0.52     |
|                                                                                                                   |
| Parallel                                                                                                          |
|   allgather()                         21        0.0865      0.004120    0.0865      0.004120    1.09     1.09     |
|   barrier()                           8         0.0001      0.000007    0.0001      0.000007    0.00     0.00     |
|   broadcast()                         25        0.0006      0.000024    0.0006      0.000024    0.01     0.01     |
|   gather()                            49        0.0118      0.000241    0.0118      0.000241    0.15     0.15     |
|   max(scalar)                         3         0.0174      0.005793    0.0174      0.005793    0.22     0.22     |
|   max(vector)                         5         0.0003      0.000069    0.0003      0.000069    0.00     0.00     |
|   maxloc(scalar)                      12        0.0807      0.006722    0.0807      0.006722    1.02     1.02     |
|   min(vector)                         5         0.0026      0.000512    0.0026      0.000512    0.03     0.03     |
|   probe()                             110       0.0412      0.000375    0.0412      0.000375    0.52     0.52     |
|   receive()                           398       0.0005      0.000001    0.0418      0.000105    0.01     0.53     |
|   send()                              158       0.0003      0.000002    0.0003      0.000002    0.00     0.00     |
|   send_receive()                      122       0.0002      0.000002    0.0422      0.000346    0.00     0.53     |
|   sum()                               24        0.0258      0.001076    0.0258      0.001076    0.33     0.33     |
|                                                                                                                   |
| Parallel::Request                                                                                                 |
|   wait()                              398       0.0003      0.000001    0.0003      0.000001    0.00     0.00     |
|                                                                                                                   |
| Partitioner                                                                                                       |
|   set_node_processor_ids()            1         0.0007      0.000717    0.0029      0.002944    0.01     0.04     |
|   set_parent_processor_ids()          1         0.0003      0.000287    0.0003      0.000287    0.00     0.00     |
|                                                                                                                   |
| PetscLinearSolver                                                                                                 |
|   solve()                             48        4.1031      0.085482    4.1031      0.085482    51.79    51.79    |
|                                                                                                                   |
| PointLocatorTree                                                                                                  |
|   init(no master)                     1         0.0059      0.005950    0.0065      0.006538    0.08     0.08     |
|   operator()                          15        0.0003      0.000021    0.0004      0.000025    0.00     0.00     |
|                                                                                                                   |
| RBConstruction                                                                                                    |
|   add_scaled_matrix_and_vector()      10        0.8418      0.084175    1.2969      0.129691    10.62    16.37    |
|   clear()                             3         0.0007      0.000240    0.0007      0.000240    0.01     0.01     |
|   compute_Fq_representor_innerprods() 2         0.1399      0.069962    0.3185      0.159262    1.77     4.02     |
|   compute_max_error_bound()           8         0.0034      0.000431    0.3203      0.040041    0.04     4.04     |
|   enrich_RB_space()                   3         0.0028      0.000934    0.0028      0.000934    0.04     0.04     |
|   train_reduced_basis()               2         0.0035      0.001775    2.6334      1.316704    0.04     33.24    |
|   truth_assembly()                    3         0.0408      0.013614    0.0409      0.013630    0.52     0.52     |
|   truth_solve()                       3         0.0023      0.000769    0.4268      0.142261    0.03     5.39     |
|   update_RB_system_matrices()         7         0.0248      0.003549    0.0248      0.003549    0.31     0.31     |
|   update_residual_terms()             3         0.0928      0.030924    1.3701      0.456695    1.17     17.29    |
|                                                                                                                   |
| RBEIMConstruction                                                                                                 |
|   compute_best_fit_error()            100       0.1690      0.001690    0.1935      0.001935    2.13     2.44     |
|   enrich_RB_space()                   4         0.0398      0.009947    0.1528      0.038203    0.50     1.93     |
|   truth_solve()                       129       0.9258      0.007177    3.5800      0.027752    11.69    45.19    |
|   update_RB_system_matrices()         4         0.0119      0.002981    0.0188      0.004708    0.15     0.24     |
|                                                                                                                   |
| RBEIMEvaluation                                                                                                   |
|   rb_solve()                          673       0.0030      0.000004    0.0030      0.000004    0.04     0.04     |
|   write_offline_data_to_files()       1         0.0002      0.000177    0.0447      0.044741    0.00     0.56     |
|                                                                                                                   |
| RBEvaluation                                                                                                      |
|   clear()                             3         0.0000      0.000015    0.0000      0.000015    0.00     0.00     |
|   compute_residual_dual_norm()        668       0.0918      0.000137    0.0918      0.000137    1.16     1.16     |
|   rb_solve()                          668       0.0055      0.000008    0.1002      0.000150    0.07     1.27     |
|   resize_data_structures()            2         0.0002      0.000084    0.0002      0.000084    0.00     0.00     |
|   write_offline_data_to_files()       2         0.0804      0.040192    0.0804      0.040192    1.01     1.01     |
 -------------------------------------------------------------------------------------------------------------------
| Totals:                               390899    7.9225                                          100.00            |
 -------------------------------------------------------------------------------------------------------------------

*** Warning, This code is untested, experimental, or likely to see future API changes: src/reduced_basis/rb_parametrized.C, line 40, compiled Aug 24 2012 at 15:15:42 ***
 EquationSystems
  n_systems()=2
   System #0, "EIM"
    Type "RBConstruction"
    Variables="x_comp_of_G" "y_comp_of_G" "z_comp_of_G" 
    Finite Element Types="LAGRANGE", "JACOBI_20_00" "LAGRANGE", "JACOBI_20_00" "LAGRANGE", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" "CARTESIAN" "CARTESIAN" 
    Approximation Orders="FIRST", "THIRD" "FIRST", "THIRD" "FIRST", "THIRD" 
    n_dofs()=18513
    n_local_dofs()=3414
    n_constrained_dofs()=0
    n_local_constrained_dofs()=0
    n_vectors()=1
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 22.0896
      Average Off-Processor Bandwidth <= 1.72056
      Maximum  On-Processor Bandwidth <= 27
      Maximum Off-Processor Bandwidth <= 14
    DofMap Constraints
      Number of DoF Constraints = 0
      Number of Node Constraints = 0
   System #1, "RB"
    Type "RBConstruction"
    Variables="u" 
    Finite Element Types="LAGRANGE", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" 
    Approximation Orders="FIRST", "THIRD" 
    n_dofs()=6171
    n_local_dofs()=1138
    n_constrained_dofs()=243
    n_local_constrained_dofs()=1
    n_vectors()=1
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 22.0896
      Average Off-Processor Bandwidth <= 1.72056
      Maximum  On-Processor Bandwidth <= 27
      Maximum Off-Processor Bandwidth <= 14
    DofMap Constraints
      Number of DoF Constraints = 242
      Average DoF Constraint Length= 0
      Number of Node Constraints = 0

 Mesh Information:
  mesh_dimension()=3
  spatial_dimension()=3
  n_nodes()=6171
    n_local_nodes()=1138
  n_elem()=5000
    n_local_elem()=833
    n_active_elem()=5000
  n_subdomains()=1
  n_partitions()=6
  n_processors()=6
  n_threads()=1
  processor_id()=0

Bi: 0.005
curvature: 1.0472
kappa: 1.3

************************************************************************************************************************
***             WIDEN YOUR WINDOW TO 120 CHARACTERS.  Use 'enscript -r -fCourier9' to print this document            ***
************************************************************************************************************************

---------------------------------------------- PETSc Performance Summary: ----------------------------------------------

./reduced_basis_ex6-opt on a intel-11. named daedalus with 6 processors, by roystgnr Fri Aug 24 15:26:10 2012
Using Petsc Release Version 3.1.0, Patch 5, Mon Sep 27 11:51:54 CDT 2010

                         Max       Max/Min        Avg      Total 
Time (sec):           5.505e-01      1.01014   5.460e-01
Objects:              3.000e+01      1.00000   3.000e+01
Flops:                2.760e+04      1.34033   2.468e+04  1.481e+05
Flops/sec:            5.064e+04      1.34063   4.520e+04  2.712e+05
MPI Messages:         1.200e+01      2.00000   1.000e+01  6.000e+01
MPI Message Lengths:  2.090e+04      2.18030   1.626e+03  9.757e+04
MPI Reductions:       6.600e+01      1.00000

Flop counting convention: 1 flop = 1 real number operation of type (multiply/divide/add/subtract)
                            e.g., VecAXPY() for real vectors of length N --> 2N flops
                            and VecAXPY() for complex vectors of length N --> 8N flops

Summary of Stages:   ----- Time ------  ----- Flops -----  --- Messages ---  -- Message Lengths --  -- Reductions --
                        Avg     %Total     Avg     %Total   counts   %Total     Avg         %Total   counts   %Total 
 0:      Main Stage: 5.4595e-01 100.0%  1.4810e+05 100.0%  6.000e+01 100.0%  1.626e+03      100.0%  3.800e+01  57.6% 

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

VecCopy                2 1.0 5.2214e-05 6.6 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecSet                16 1.0 5.8413e-05 2.3 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecAXPY                6 1.0 1.0490e-04 1.6 2.76e+04 1.3 0.0e+00 0.0e+00 0.0e+00  0100  0  0  0   0100  0  0  0  1412
VecAssemblyBegin       8 1.0 4.6558e-03 1.7 0.00e+00 0.0 0.0e+00 0.0e+00 2.4e+01  1  0  0  0 36   1  0  0  0 63     0
VecAssemblyEnd         8 1.0 5.5075e-05 1.4 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecScatterBegin        2 1.0 7.3910e-05 2.2 0.00e+00 0.0 2.0e+01 3.2e+03 0.0e+00  0  0 33 67  0   0  0 33 67  0     0
VecScatterEnd          2 1.0 1.5283e-0413.9 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatZeroEntries         4 1.0 4.6682e-04 2.3 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
------------------------------------------------------------------------------------------------------------------------

Memory usage is given in bytes:

Object Type          Creations   Destructions     Memory  Descendants' Mem.
Reports information only for process 0.

--- Event Stage 0: Main Stage

                 Vec    16             16       284480     0
         Vec Scatter     2              2         1736     0
           Index Set     4              4         6160     0
   IS L to G Mapping     2              2        23112     0
              Matrix     6              6      1422928     0
========================================================================================================================
Average time to get PetscTime(): 1.19209e-07
Average time for MPI_Barrier(): 4.69685e-05
Average time for zero size MPI_Send(): 4.38293e-05
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
| Time:           Fri Aug 24 15:26:10 2012                         |
| OS:             Linux                                            |
| HostName:       daedalus                                         |
| OS Release:     2.6.32-34-generic                                |
| OS Version:     #76-Ubuntu SMP Tue Aug 30 17:05:01 UTC 2011      |
| Machine:        x86_64                                           |
| Username:       roystgnr                                         |
| Configuration:  ./configure run on Wed Aug 22 12:44:06 CDT 2012  |
-------------------------------------------------------------------
 --------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=0.629518, Active time=0.457937                                               |
 --------------------------------------------------------------------------------------------------------------
| Event                            nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                            w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|--------------------------------------------------------------------------------------------------------------|
|                                                                                                              |
|                                                                                                              |
| DofMap                                                                                                       |
|   add_neighbors_to_send_list()   2         0.0016      0.000778    0.0020      0.001010    0.34     0.44     |
|   build_sparsity()               2         0.0279      0.013939    0.0361      0.018038    6.09     7.88     |
|   create_dof_constraints()       2         0.0086      0.004291    0.0106      0.005306    1.87     2.32     |
|   distribute_dofs()              2         0.0045      0.002229    0.0382      0.019083    0.97     8.33     |
|   dof_indices()                  29184     0.0102      0.000000    0.0102      0.000000    2.22     2.22     |
|   prepare_send_list()            2         0.0002      0.000121    0.0002      0.000121    0.05     0.05     |
|   reinit()                       2         0.0130      0.006517    0.0130      0.006517    2.85     2.85     |
|                                                                                                              |
| EquationSystems                                                                                              |
|   build_solution_vector()        1         0.0037      0.003690    0.0232      0.023194    0.81     5.06     |
|                                                                                                              |
| ExodusII_IO                                                                                                  |
|   write_nodal_data()             1         0.0063      0.006299    0.0063      0.006299    1.38     1.38     |
|                                                                                                              |
| Mesh                                                                                                         |
|   find_neighbors()               1         0.0083      0.008253    0.0198      0.019795    1.80     4.32     |
|   renumber_nodes_and_elem()      2         0.0007      0.000374    0.0007      0.000374    0.16     0.16     |
|                                                                                                              |
| MeshCommunication                                                                                            |
|   assign_global_indices()        2         0.1011      0.050573    0.1993      0.099669    22.09    43.53    |
|   compute_hilbert_indices()      2         0.0193      0.009629    0.0193      0.009629    4.21     4.21     |
|   find_global_indices()          2         0.0020      0.000999    0.0459      0.022970    0.44     10.03    |
|   parallel_sort()                2         0.0017      0.000842    0.0185      0.009257    0.37     4.04     |
|                                                                                                              |
| MeshOutput                                                                                                   |
|   write_equation_systems()       1         0.0000      0.000031    0.0295      0.029524    0.01     6.45     |
|                                                                                                              |
| MeshTools::Generation                                                                                        |
|   build_cube()                   1         0.0022      0.002217    0.0022      0.002217    0.48     0.48     |
|                                                                                                              |
| MetisPartitioner                                                                                             |
|   partition()                    1         0.0103      0.010335    0.0354      0.035430    2.26     7.74     |
|                                                                                                              |
| Parallel                                                                                                     |
|   allgather()                    45        0.0951      0.002114    0.0951      0.002114    20.78    20.78    |
|   barrier()                      2         0.0000      0.000012    0.0000      0.000012    0.01     0.01     |
|   broadcast()                    223       0.0015      0.000007    0.0015      0.000007    0.32     0.32     |
|   gather()                       1         0.0000      0.000004    0.0000      0.000004    0.00     0.00     |
|   max(scalar)                    3         0.0243      0.008101    0.0243      0.008101    5.31     5.31     |
|   max(vector)                    4         0.0002      0.000044    0.0002      0.000044    0.04     0.04     |
|   min(vector)                    4         0.0027      0.000670    0.0027      0.000670    0.59     0.59     |
|   probe()                        110       0.0545      0.000496    0.0545      0.000496    11.91    11.91    |
|   receive()                      110       0.0004      0.000004    0.0550      0.000500    0.09     12.00    |
|   send()                         110       0.0003      0.000002    0.0003      0.000002    0.06     0.06     |
|   send_receive()                 122       0.0003      0.000002    0.0556      0.000456    0.06     12.15    |
|   sum()                          24        0.0436      0.001817    0.0436      0.001817    9.52     9.52     |
|                                                                                                              |
| Parallel::Request                                                                                            |
|   wait()                         110       0.0001      0.000001    0.0001      0.000001    0.01     0.01     |
|                                                                                                              |
| Partitioner                                                                                                  |
|   set_node_processor_ids()       1         0.0007      0.000710    0.0105      0.010502    0.16     2.29     |
|   set_parent_processor_ids()     1         0.0003      0.000287    0.0003      0.000287    0.06     0.06     |
|                                                                                                              |
| RBConstruction                                                                                               |
|   clear()                        3         0.0004      0.000128    0.0004      0.000128    0.08     0.08     |
|   load_rb_solution()             2         0.0009      0.000425    0.0009      0.000425    0.19     0.19     |
|                                                                                                              |
| RBEIMEvaluation                                                                                              |
|   rb_solve()                     1         0.0100      0.009959    0.0100      0.009959    2.17     2.17     |
|   read_offline_data_from_files() 1         0.0001      0.000055    0.0006      0.000644    0.01     0.14     |
|                                                                                                              |
| RBEvaluation                                                                                                 |
|   clear()                        3         0.0001      0.000025    0.0001      0.000025    0.02     0.02     |
|   compute_residual_dual_norm()   1         0.0002      0.000204    0.0002      0.000204    0.04     0.04     |
|   rb_solve()                     1         0.0000      0.000032    0.0102      0.010195    0.01     2.23     |
|   read_offline_data_from_files() 2         0.0007      0.000369    0.0008      0.000400    0.16     0.17     |
|   resize_data_structures()       2         0.0001      0.000030    0.0001      0.000030    0.01     0.01     |
 --------------------------------------------------------------------------------------------------------------
| Totals:                          30098     0.4579                                          100.00            |
 --------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running  mpirun -np 6 ./reduced_basis_ex6-opt -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
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
