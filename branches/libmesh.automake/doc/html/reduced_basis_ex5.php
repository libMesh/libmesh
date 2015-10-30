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
<h1>Reduced Basis: Example 5 - Reduced Basis Cantilever</h1>
3D cantilever beam using the Reduced Basis Method
David Knezevic, Kyunghoon "K" Lee

<br><br>We consider four parameters in this problem:
x_scaling: scales the length of the cantilever
load_Fx: the traction in the x-direction on the right boundary of the cantilever
load_Fy: the traction in the y-direction on the right boundary of the cantilever
load_Fz: the traction in the z-direction on the right boundary of the cantilever


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
        #include "quadrature_gauss.h"
        #include "libmesh_logging.h"
        
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
Define a function to scale the mesh according to the parameter.
</div>

<div class ="fragment">
<pre>
        void scale_mesh_and_plot(EquationSystems& es, const RBParameters& mu, const std::string& filename);
        
</pre>
</div>
<div class = "comment">
Post-process the solution to compute stresses
</div>

<div class ="fragment">
<pre>
        void compute_stresses(EquationSystems& es);
        
</pre>
</div>
<div class = "comment">
The main program.
</div>

<div class ="fragment">
<pre>
        int main(int argc, char** argv) {
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
This example only works if libMesh was compiled for 3D
</div>

<div class ="fragment">
<pre>
          const unsigned int dim = 3;
          libmesh_example_assert(dim == LIBMESH_DIM, "3D support");
        
        	std::string parameters_filename = "reduced_basis_ex5.in";
        	GetPot infile(parameters_filename);
        
          unsigned int n_elem_x  = infile("n_elem_x",0);
          unsigned int n_elem_y  = infile("n_elem_y",0);
          unsigned int n_elem_z  = infile("n_elem_z",0);
          Real x_size            = infile("x_size", 0.);
          Real y_size            = infile("y_size", 0.);
          Real z_size            = infile("z_size", 0.);
        
        	bool store_basis_functions = infile("store_basis_functions", true);
        
</pre>
</div>
<div class = "comment">
Read the "online_mode" flag from the command line
</div>

<div class ="fragment">
<pre>
                GetPot command_line(argc, argv);
        	int online_mode = 0;
        	if ( command_line.search(1, "-online_mode") ) {
        		online_mode = command_line.next(online_mode);		
        	}
        
        
          Mesh mesh (dim);
          MeshTools::Generation::build_cube (mesh,
                                             n_elem_x,
                                             n_elem_y,
                                             n_elem_z,
                                             0., x_size,
                                             0., y_size,
                                             0., z_size,
                                             HEX8);
        	 mesh.print_info();
        
</pre>
</div>
<div class = "comment">
Create an equation systems object.
</div>

<div class ="fragment">
<pre>
                EquationSystems equation_systems(mesh);
        
</pre>
</div>
<div class = "comment">
We override RBConstruction with ElasticityRBConstruction in order to
specialize a few functions for this particular problem.
</div>

<div class ="fragment">
<pre>
                ElasticityRBConstruction& rb_con =
        		equation_systems.add_system&lt;ElasticityRBConstruction&gt;("RBElasticity");
        
</pre>
</div>
<div class = "comment">
Also, initialize an ExplicitSystem to store stresses
</div>

<div class ="fragment">
<pre>
          ExplicitSystem& stress_system =
            equation_systems.add_system&lt;ExplicitSystem&gt; ("StressSystem");
          stress_system.add_variable("sigma_00", CONSTANT, MONOMIAL);
          stress_system.add_variable("sigma_01", CONSTANT, MONOMIAL);
          stress_system.add_variable("sigma_02", CONSTANT, MONOMIAL);
          stress_system.add_variable("sigma_10", CONSTANT, MONOMIAL);
          stress_system.add_variable("sigma_11", CONSTANT, MONOMIAL);
          stress_system.add_variable("sigma_12", CONSTANT, MONOMIAL);
          stress_system.add_variable("sigma_20", CONSTANT, MONOMIAL);
          stress_system.add_variable("sigma_21", CONSTANT, MONOMIAL);
          stress_system.add_variable("sigma_22", CONSTANT, MONOMIAL);
          stress_system.add_variable("vonMises", CONSTANT, MONOMIAL);
        
</pre>
</div>
<div class = "comment">
Initialize the data structures for the equation system.
</div>

<div class ="fragment">
<pre>
                equation_systems.init ();
          equation_systems.print_info();
        
</pre>
</div>
<div class = "comment">
Build a new RBEvaluation object which will be used to perform
Reduced Basis calculations. This is required in both the
"Offline" and "Online" stages.
</div>

<div class ="fragment">
<pre>
                ElasticityRBEvaluation rb_eval;
        
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
Iinitialize online parameters
</div>

<div class ="fragment">
<pre>
                        Real online_x_scaling = infile("online_x_scaling", 0.);
        		Real online_load_Fx   = infile("online_load_Fx",   0.);
        		Real online_load_Fy   = infile("online_load_Fy",   0.);
        		Real online_load_Fz   = infile("online_load_Fz",   0.);
        		RBParameters online_mu;
        		online_mu.set_value("x_scaling", online_x_scaling);
        		online_mu.set_value("load_Fx",   online_load_Fx);
        		online_mu.set_value("load_Fy",   online_load_Fy);
        		online_mu.set_value("load_Fz",   online_load_Fz);
        		rb_eval.set_parameters(online_mu);
        		rb_eval.print_parameters();
        		
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
        
        			const RBParameters& rb_eval_params = rb_eval.get_parameters();
        			scale_mesh_and_plot(equation_systems, rb_eval_params, "RB_sol.e");
        		}
        	}
        
        	return 0;
        }
        
        void scale_mesh_and_plot(EquationSystems& es, const RBParameters& mu, const std::string& filename)
        {
          const Real x_scaling = mu.get_value("x_scaling");
          
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
        
            (*node)(0) *= mu.get_value("x_scaling");
          }
          
</pre>
</div>
<div class = "comment">
Post-process the solution to compute the stresses
</div>

<div class ="fragment">
<pre>
          compute_stresses(es);
        
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
          node_it = mesh.nodes_begin();
        
          for( ; node_it != node_end; node_it++)
          {
            Node* node = *node_it;
        
            (*node)(0) /= mu.get_value("x_scaling");
          }
        }
        
        void compute_stresses(EquationSystems& es)
        {
          START_LOG("compute_stresses()", "main");
        
          const MeshBase& mesh = es.get_mesh();
        
          const unsigned int dim = mesh.mesh_dimension();
        
          ElasticityRBConstruction& system = es.get_system&lt;ElasticityRBConstruction&gt;("RBElasticity");
        
          unsigned int displacement_vars[3];
          displacement_vars[0] = system.variable_number ("u");
          displacement_vars[1] = system.variable_number ("v");
          displacement_vars[2] = system.variable_number ("w");
          const unsigned int u_var = system.variable_number ("u");
        
          const DofMap& dof_map = system.get_dof_map();
          FEType fe_type = dof_map.variable_type(u_var);
          AutoPtr&lt;FEBase&gt; fe (FEBase::build(dim, fe_type));
          QGauss qrule (dim, fe_type.default_quadrature_order());
          fe-&gt;attach_quadrature_rule (&qrule);
          
          const std::vector&lt;Real&gt;& JxW = fe-&gt;get_JxW();
          const std::vector&lt;std::vector&lt;RealGradient&gt; &gt;& dphi = fe-&gt;get_dphi();
          
</pre>
</div>
<div class = "comment">
Also, get a reference to the ExplicitSystem
</div>

<div class ="fragment">
<pre>
          ExplicitSystem& stress_system = es.get_system&lt;ExplicitSystem&gt;("StressSystem");
          const DofMap& stress_dof_map = stress_system.get_dof_map();
          unsigned int sigma_vars[3][3];
          sigma_vars[0][0] = stress_system.variable_number ("sigma_00");
          sigma_vars[0][1] = stress_system.variable_number ("sigma_01");
          sigma_vars[0][2] = stress_system.variable_number ("sigma_02");
          sigma_vars[1][0] = stress_system.variable_number ("sigma_10");
          sigma_vars[1][1] = stress_system.variable_number ("sigma_11");
          sigma_vars[1][2] = stress_system.variable_number ("sigma_12");
          sigma_vars[2][0] = stress_system.variable_number ("sigma_20");
          sigma_vars[2][1] = stress_system.variable_number ("sigma_21");
          sigma_vars[2][2] = stress_system.variable_number ("sigma_22");
          unsigned int vonMises_var = stress_system.variable_number ("vonMises");
        
</pre>
</div>
<div class = "comment">
Storage for the stress dof indices on each element
</div>

<div class ="fragment">
<pre>
          std::vector&lt; std::vector&lt;unsigned int&gt; &gt; dof_indices_var(system.n_vars());
          std::vector&lt;unsigned int&gt; stress_dof_indices_var;
        
</pre>
</div>
<div class = "comment">
To store the stress tensor on each element
</div>

<div class ="fragment">
<pre>
          DenseMatrix&lt;Number&gt; elem_sigma;
        
          MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
          const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
        
          for ( ; el != end_el; ++el)
          {
            const Elem* elem = *el;
        
            for(unsigned int var=0; var&lt;3; var++)
            {
              dof_map.dof_indices (elem, dof_indices_var[var], displacement_vars[var]);
            }
        
            fe-&gt;reinit (elem);
        
            elem_sigma.resize(3,3);
            
            for (unsigned int qp=0; qp&lt;qrule.n_points(); qp++)
            {
              for(unsigned int C_i=0; C_i&lt;3; C_i++)
                for(unsigned int C_j=0; C_j&lt;3; C_j++)
                  for(unsigned int C_k=0; C_k&lt;3; C_k++)
                  {
                    const unsigned int n_var_dofs = dof_indices_var[C_k].size();
        
</pre>
</div>
<div class = "comment">
Get the gradient at this quadrature point
</div>

<div class ="fragment">
<pre>
                    Gradient displacement_gradient;
                    for(unsigned int l=0; l&lt;n_var_dofs; l++)
                    {
                      displacement_gradient.add_scaled(dphi[l][qp], system.current_solution(dof_indices_var[C_k][l]));
                    }
        
                    for(unsigned int C_l=0; C_l&lt;3; C_l++)
                    {
                      elem_sigma(C_i,C_j) += JxW[qp]*( elasticity_tensor(C_i,C_j,C_k,C_l) * displacement_gradient(C_l) );
                    }
        
                  }
            }
            
</pre>
</div>
<div class = "comment">
Get the average stresses by dividing by the element volume
</div>

<div class ="fragment">
<pre>
            elem_sigma.scale(1./elem-&gt;volume());
        
</pre>
</div>
<div class = "comment">
load elem_sigma data into stress_system
</div>

<div class ="fragment">
<pre>
            for(unsigned int i=0; i&lt;3; i++)
              for(unsigned int j=0; j&lt;3; j++)
              {
                stress_dof_map.dof_indices (elem, stress_dof_indices_var, sigma_vars[i][j]);
        
</pre>
</div>
<div class = "comment">
We are using CONSTANT MONOMIAL basis functions, hence we only need to get
one dof index per variable
</div>

<div class ="fragment">
<pre>
                unsigned int dof_index = stress_dof_indices_var[0];
                
                if( (stress_system.solution-&gt;first_local_index() &lt;= dof_index) &&
                    (dof_index &lt; stress_system.solution-&gt;last_local_index()) )
                {
                  stress_system.solution-&gt;set(dof_index, elem_sigma(i,j));
                }
        
              }
            
</pre>
</div>
<div class = "comment">
Also, the von Mises stress
</div>

<div class ="fragment">
<pre>
            Number vonMises_value = std::sqrt( 0.5*( pow(elem_sigma(0,0) - elem_sigma(1,1),2.) + 
                                                     pow(elem_sigma(1,1) - elem_sigma(2,2),2.) + 
                                                     pow(elem_sigma(2,2) - elem_sigma(0,0),2.) +
                                                     6.*(pow(elem_sigma(0,1),2.) + pow(elem_sigma(1,2),2.) + pow(elem_sigma(2,0),2.))
                                                   ) );
            stress_dof_map.dof_indices (elem, stress_dof_indices_var, vonMises_var);
            unsigned int dof_index = stress_dof_indices_var[0];
            if( (stress_system.solution-&gt;first_local_index() &lt;= dof_index) &&
                (dof_index &lt; stress_system.solution-&gt;last_local_index()) )
            {
              stress_system.solution-&gt;set(dof_index, vonMises_value);
            }
            
          }
        
</pre>
</div>
<div class = "comment">
Should call close and update when we set vector entries directly
</div>

<div class ="fragment">
<pre>
          stress_system.solution-&gt;close();
          stress_system.update();
        
          STOP_LOG("compute_stresses()", "main");
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
  #include <B><FONT COLOR="#BC8F8F">&quot;quadrature_gauss.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh_logging.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;rb_classes.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;assembly.h&quot;</FONT></B>
  
  using namespace libMesh;
  
  <B><FONT COLOR="#228B22">void</FONT></B> scale_mesh_and_plot(EquationSystems&amp; es, <B><FONT COLOR="#228B22">const</FONT></B> RBParameters&amp; mu, <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp; filename);
  
  <B><FONT COLOR="#228B22">void</FONT></B> compute_stresses(EquationSystems&amp; es);
  
  <B><FONT COLOR="#228B22">int</FONT></B> main(<B><FONT COLOR="#228B22">int</FONT></B> argc, <B><FONT COLOR="#228B22">char</FONT></B>** argv) {
  	LibMeshInit init (argc, argv);
  
  #<B><FONT COLOR="#A020F0">if</FONT></B> !defined(LIBMESH_HAVE_XDR)
    libmesh_example_assert(false, <B><FONT COLOR="#BC8F8F">&quot;--enable-xdr&quot;</FONT></B>);
  #elif defined(LIBMESH_DEFAULT_SINGLE_PRECISION)
    libmesh_example_assert(false, <B><FONT COLOR="#BC8F8F">&quot;--disable-singleprecision&quot;</FONT></B>);
  #endif
  
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> dim = 3;
    libmesh_example_assert(dim == LIBMESH_DIM, <B><FONT COLOR="#BC8F8F">&quot;3D support&quot;</FONT></B>);
  
  	<B><FONT COLOR="#5F9EA0">std</FONT></B>::string parameters_filename = <B><FONT COLOR="#BC8F8F">&quot;reduced_basis_ex5.in&quot;</FONT></B>;
  	GetPot infile(parameters_filename);
  
    <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_elem_x  = infile(<B><FONT COLOR="#BC8F8F">&quot;n_elem_x&quot;</FONT></B>,0);
    <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_elem_y  = infile(<B><FONT COLOR="#BC8F8F">&quot;n_elem_y&quot;</FONT></B>,0);
    <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_elem_z  = infile(<B><FONT COLOR="#BC8F8F">&quot;n_elem_z&quot;</FONT></B>,0);
    Real x_size            = infile(<B><FONT COLOR="#BC8F8F">&quot;x_size&quot;</FONT></B>, 0.);
    Real y_size            = infile(<B><FONT COLOR="#BC8F8F">&quot;y_size&quot;</FONT></B>, 0.);
    Real z_size            = infile(<B><FONT COLOR="#BC8F8F">&quot;z_size&quot;</FONT></B>, 0.);
  
  	<B><FONT COLOR="#228B22">bool</FONT></B> store_basis_functions = infile(<B><FONT COLOR="#BC8F8F">&quot;store_basis_functions&quot;</FONT></B>, true);
  
  	GetPot command_line(argc, argv);
  	<B><FONT COLOR="#228B22">int</FONT></B> online_mode = 0;
  	<B><FONT COLOR="#A020F0">if</FONT></B> ( command_line.search(1, <B><FONT COLOR="#BC8F8F">&quot;-online_mode&quot;</FONT></B>) ) {
  		online_mode = command_line.next(online_mode);		
  	}
  
  
    Mesh mesh (dim);
    <B><FONT COLOR="#5F9EA0">MeshTools</FONT></B>::Generation::build_cube (mesh,
                                       n_elem_x,
                                       n_elem_y,
                                       n_elem_z,
                                       0., x_size,
                                       0., y_size,
                                       0., z_size,
                                       HEX8);
  	 mesh.print_info();
  
  	EquationSystems equation_systems(mesh);
  
  	ElasticityRBConstruction&amp; rb_con =
  		equation_systems.add_system&lt;ElasticityRBConstruction&gt;(<B><FONT COLOR="#BC8F8F">&quot;RBElasticity&quot;</FONT></B>);
  
    ExplicitSystem&amp; stress_system =
      equation_systems.add_system&lt;ExplicitSystem&gt; (<B><FONT COLOR="#BC8F8F">&quot;StressSystem&quot;</FONT></B>);
    stress_system.add_variable(<B><FONT COLOR="#BC8F8F">&quot;sigma_00&quot;</FONT></B>, CONSTANT, MONOMIAL);
    stress_system.add_variable(<B><FONT COLOR="#BC8F8F">&quot;sigma_01&quot;</FONT></B>, CONSTANT, MONOMIAL);
    stress_system.add_variable(<B><FONT COLOR="#BC8F8F">&quot;sigma_02&quot;</FONT></B>, CONSTANT, MONOMIAL);
    stress_system.add_variable(<B><FONT COLOR="#BC8F8F">&quot;sigma_10&quot;</FONT></B>, CONSTANT, MONOMIAL);
    stress_system.add_variable(<B><FONT COLOR="#BC8F8F">&quot;sigma_11&quot;</FONT></B>, CONSTANT, MONOMIAL);
    stress_system.add_variable(<B><FONT COLOR="#BC8F8F">&quot;sigma_12&quot;</FONT></B>, CONSTANT, MONOMIAL);
    stress_system.add_variable(<B><FONT COLOR="#BC8F8F">&quot;sigma_20&quot;</FONT></B>, CONSTANT, MONOMIAL);
    stress_system.add_variable(<B><FONT COLOR="#BC8F8F">&quot;sigma_21&quot;</FONT></B>, CONSTANT, MONOMIAL);
    stress_system.add_variable(<B><FONT COLOR="#BC8F8F">&quot;sigma_22&quot;</FONT></B>, CONSTANT, MONOMIAL);
    stress_system.add_variable(<B><FONT COLOR="#BC8F8F">&quot;vonMises&quot;</FONT></B>, CONSTANT, MONOMIAL);
  
  	equation_systems.init ();
    equation_systems.print_info();
  
  	ElasticityRBEvaluation rb_eval;
  
  	rb_con.set_rb_evaluation(rb_eval);
  
  	<B><FONT COLOR="#A020F0">if</FONT></B>(!online_mode) <I><FONT COLOR="#B22222">// Perform the Offline stage of the RB method
</FONT></I>  	{
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
</FONT></I>  	{
  		rb_eval.read_offline_data_from_files();
  
  		Real online_x_scaling = infile(<B><FONT COLOR="#BC8F8F">&quot;online_x_scaling&quot;</FONT></B>, 0.);
  		Real online_load_Fx   = infile(<B><FONT COLOR="#BC8F8F">&quot;online_load_Fx&quot;</FONT></B>,   0.);
  		Real online_load_Fy   = infile(<B><FONT COLOR="#BC8F8F">&quot;online_load_Fy&quot;</FONT></B>,   0.);
  		Real online_load_Fz   = infile(<B><FONT COLOR="#BC8F8F">&quot;online_load_Fz&quot;</FONT></B>,   0.);
  		RBParameters online_mu;
  		online_mu.set_value(<B><FONT COLOR="#BC8F8F">&quot;x_scaling&quot;</FONT></B>, online_x_scaling);
  		online_mu.set_value(<B><FONT COLOR="#BC8F8F">&quot;load_Fx&quot;</FONT></B>,   online_load_Fx);
  		online_mu.set_value(<B><FONT COLOR="#BC8F8F">&quot;load_Fy&quot;</FONT></B>,   online_load_Fy);
  		online_mu.set_value(<B><FONT COLOR="#BC8F8F">&quot;load_Fz&quot;</FONT></B>,   online_load_Fz);
  		rb_eval.set_parameters(online_mu);
  		rb_eval.print_parameters();
  		
  		rb_eval.rb_solve( rb_eval.get_n_basis_functions() );
  
  		<B><FONT COLOR="#A020F0">if</FONT></B>(store_basis_functions)
  		{
  			rb_eval.read_in_basis_functions(rb_con);
  
  			rb_con.load_rb_solution();
  
  			<B><FONT COLOR="#228B22">const</FONT></B> RBParameters&amp; rb_eval_params = rb_eval.get_parameters();
  			scale_mesh_and_plot(equation_systems, rb_eval_params, <B><FONT COLOR="#BC8F8F">&quot;RB_sol.e&quot;</FONT></B>);
  		}
  	}
  
  	<B><FONT COLOR="#A020F0">return</FONT></B> 0;
  }
  
  <B><FONT COLOR="#228B22">void</FONT></B> scale_mesh_and_plot(EquationSystems&amp; es, <B><FONT COLOR="#228B22">const</FONT></B> RBParameters&amp; mu, <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp; filename)
  {
    <B><FONT COLOR="#228B22">const</FONT></B> Real x_scaling = mu.get_value(<B><FONT COLOR="#BC8F8F">&quot;x_scaling&quot;</FONT></B>);
    
    MeshBase&amp; mesh = es.get_mesh();
  
    <B><FONT COLOR="#5F9EA0">MeshBase</FONT></B>::node_iterator       node_it  = mesh.nodes_begin();
    <B><FONT COLOR="#228B22">const</FONT></B> MeshBase::node_iterator node_end = mesh.nodes_end();
  
    <B><FONT COLOR="#A020F0">for</FONT></B>( ; node_it != node_end; node_it++)
    {
      Node* node = *node_it;
  
      (*node)(0) *= mu.get_value(<B><FONT COLOR="#BC8F8F">&quot;x_scaling&quot;</FONT></B>);
    }
    
    compute_stresses(es);
  
  #ifdef LIBMESH_HAVE_EXODUS_API
    ExodusII_IO (mesh).write_equation_systems (filename, es);
  #endif
  
    node_it = mesh.nodes_begin();
  
    <B><FONT COLOR="#A020F0">for</FONT></B>( ; node_it != node_end; node_it++)
    {
      Node* node = *node_it;
  
      (*node)(0) /= mu.get_value(<B><FONT COLOR="#BC8F8F">&quot;x_scaling&quot;</FONT></B>);
    }
  }
  
  <B><FONT COLOR="#228B22">void</FONT></B> compute_stresses(EquationSystems&amp; es)
  {
    START_LOG(<B><FONT COLOR="#BC8F8F">&quot;compute_stresses()&quot;</FONT></B>, <B><FONT COLOR="#BC8F8F">&quot;main&quot;</FONT></B>);
  
    <B><FONT COLOR="#228B22">const</FONT></B> MeshBase&amp; mesh = es.get_mesh();
  
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> dim = mesh.mesh_dimension();
  
    ElasticityRBConstruction&amp; system = es.get_system&lt;ElasticityRBConstruction&gt;(<B><FONT COLOR="#BC8F8F">&quot;RBElasticity&quot;</FONT></B>);
  
    <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> displacement_vars[3];
    displacement_vars[0] = system.variable_number (<B><FONT COLOR="#BC8F8F">&quot;u&quot;</FONT></B>);
    displacement_vars[1] = system.variable_number (<B><FONT COLOR="#BC8F8F">&quot;v&quot;</FONT></B>);
    displacement_vars[2] = system.variable_number (<B><FONT COLOR="#BC8F8F">&quot;w&quot;</FONT></B>);
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> u_var = system.variable_number (<B><FONT COLOR="#BC8F8F">&quot;u&quot;</FONT></B>);
  
    <B><FONT COLOR="#228B22">const</FONT></B> DofMap&amp; dof_map = system.get_dof_map();
    FEType fe_type = dof_map.variable_type(u_var);
    AutoPtr&lt;FEBase&gt; fe (FEBase::build(dim, fe_type));
    QGauss qrule (dim, fe_type.default_quadrature_order());
    fe-&gt;attach_quadrature_rule (&amp;qrule);
    
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Real&gt;&amp; JxW = fe-&gt;get_JxW();
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;RealGradient&gt; &gt;&amp; dphi = fe-&gt;get_dphi();
    
    ExplicitSystem&amp; stress_system = es.get_system&lt;ExplicitSystem&gt;(<B><FONT COLOR="#BC8F8F">&quot;StressSystem&quot;</FONT></B>);
    <B><FONT COLOR="#228B22">const</FONT></B> DofMap&amp; stress_dof_map = stress_system.get_dof_map();
    <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> sigma_vars[3][3];
    sigma_vars[0][0] = stress_system.variable_number (<B><FONT COLOR="#BC8F8F">&quot;sigma_00&quot;</FONT></B>);
    sigma_vars[0][1] = stress_system.variable_number (<B><FONT COLOR="#BC8F8F">&quot;sigma_01&quot;</FONT></B>);
    sigma_vars[0][2] = stress_system.variable_number (<B><FONT COLOR="#BC8F8F">&quot;sigma_02&quot;</FONT></B>);
    sigma_vars[1][0] = stress_system.variable_number (<B><FONT COLOR="#BC8F8F">&quot;sigma_10&quot;</FONT></B>);
    sigma_vars[1][1] = stress_system.variable_number (<B><FONT COLOR="#BC8F8F">&quot;sigma_11&quot;</FONT></B>);
    sigma_vars[1][2] = stress_system.variable_number (<B><FONT COLOR="#BC8F8F">&quot;sigma_12&quot;</FONT></B>);
    sigma_vars[2][0] = stress_system.variable_number (<B><FONT COLOR="#BC8F8F">&quot;sigma_20&quot;</FONT></B>);
    sigma_vars[2][1] = stress_system.variable_number (<B><FONT COLOR="#BC8F8F">&quot;sigma_21&quot;</FONT></B>);
    sigma_vars[2][2] = stress_system.variable_number (<B><FONT COLOR="#BC8F8F">&quot;sigma_22&quot;</FONT></B>);
    <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> vonMises_var = stress_system.variable_number (<B><FONT COLOR="#BC8F8F">&quot;vonMises&quot;</FONT></B>);
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt; std::vector&lt;<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B>&gt; &gt; dof_indices_var(system.n_vars());
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B>&gt; stress_dof_indices_var;
  
    DenseMatrix&lt;Number&gt; elem_sigma;
  
    <B><FONT COLOR="#5F9EA0">MeshBase</FONT></B>::const_element_iterator       el     = mesh.active_local_elements_begin();
    <B><FONT COLOR="#228B22">const</FONT></B> MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
  
    <B><FONT COLOR="#A020F0">for</FONT></B> ( ; el != end_el; ++el)
    {
      <B><FONT COLOR="#228B22">const</FONT></B> Elem* elem = *el;
  
      <B><FONT COLOR="#A020F0">for</FONT></B>(<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> var=0; var&lt;3; var++)
      {
        dof_map.dof_indices (elem, dof_indices_var[var], displacement_vars[var]);
      }
  
      fe-&gt;reinit (elem);
  
      elem_sigma.resize(3,3);
      
      <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> qp=0; qp&lt;qrule.n_points(); qp++)
      {
        <B><FONT COLOR="#A020F0">for</FONT></B>(<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> C_i=0; C_i&lt;3; C_i++)
          <B><FONT COLOR="#A020F0">for</FONT></B>(<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> C_j=0; C_j&lt;3; C_j++)
            <B><FONT COLOR="#A020F0">for</FONT></B>(<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> C_k=0; C_k&lt;3; C_k++)
            {
              <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_var_dofs = dof_indices_var[C_k].size();
  
              Gradient displacement_gradient;
              <B><FONT COLOR="#A020F0">for</FONT></B>(<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> l=0; l&lt;n_var_dofs; l++)
              {
                displacement_gradient.add_scaled(dphi[l][qp], system.current_solution(dof_indices_var[C_k][l]));
              }
  
              <B><FONT COLOR="#A020F0">for</FONT></B>(<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> C_l=0; C_l&lt;3; C_l++)
              {
                elem_sigma(C_i,C_j) += JxW[qp]*( elasticity_tensor(C_i,C_j,C_k,C_l) * displacement_gradient(C_l) );
              }
  
            }
      }
      
      elem_sigma.scale(1./elem-&gt;volume());
  
      <B><FONT COLOR="#A020F0">for</FONT></B>(<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i&lt;3; i++)
        <B><FONT COLOR="#A020F0">for</FONT></B>(<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j=0; j&lt;3; j++)
        {
          stress_dof_map.dof_indices (elem, stress_dof_indices_var, sigma_vars[i][j]);
  
          <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> dof_index = stress_dof_indices_var[0];
          
          <B><FONT COLOR="#A020F0">if</FONT></B>( (stress_system.solution-&gt;first_local_index() &lt;= dof_index) &amp;&amp;
              (dof_index &lt; stress_system.solution-&gt;last_local_index()) )
          {
            stress_system.solution-&gt;set(dof_index, elem_sigma(i,j));
          }
  
        }
      
      Number vonMises_value = std::sqrt( 0.5*( pow(elem_sigma(0,0) - elem_sigma(1,1),2.) + 
                                               pow(elem_sigma(1,1) - elem_sigma(2,2),2.) + 
                                               pow(elem_sigma(2,2) - elem_sigma(0,0),2.) +
                                               6.*(pow(elem_sigma(0,1),2.) + pow(elem_sigma(1,2),2.) + pow(elem_sigma(2,0),2.))
                                             ) );
      stress_dof_map.dof_indices (elem, stress_dof_indices_var, vonMises_var);
      <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> dof_index = stress_dof_indices_var[0];
      <B><FONT COLOR="#A020F0">if</FONT></B>( (stress_system.solution-&gt;first_local_index() &lt;= dof_index) &amp;&amp;
          (dof_index &lt; stress_system.solution-&gt;last_local_index()) )
      {
        stress_system.solution-&gt;set(dof_index, vonMises_value);
      }
      
    }
  
    stress_system.solution-&gt;close();
    stress_system.update();
  
    STOP_LOG(<B><FONT COLOR="#BC8F8F">&quot;compute_stresses()&quot;</FONT></B>, <B><FONT COLOR="#BC8F8F">&quot;main&quot;</FONT></B>);
  }
</pre> 
<a name="output"></a> 
<br><br><br> <h1> The console output of the program: </h1> 
<pre>
Linking reduced_basis_ex5-opt...
***************************************************************
* Running  mpirun -np 6 ./reduced_basis_ex5-opt -ksp_type cg -online_mode ?
***************************************************************
 
 Mesh Information:
  mesh_dimension()=3
  spatial_dimension()=3
  n_nodes()=1845
    n_local_nodes()=348
  n_elem()=1280
    n_local_elem()=214
    n_active_elem()=1280
  n_subdomains()=1
  n_partitions()=6
  n_processors()=6
  n_threads()=1
  processor_id()=0

*** Warning, This code is untested, experimental, or likely to see future API changes: src/reduced_basis/rb_parametrized.C, line 40, compiled Aug 24 2012 at 15:15:42 ***
 EquationSystems
  n_systems()=2
   System #0, "RBElasticity"
    Type "RBConstruction"
    Variables="u" "v" "w" 
    Finite Element Types="LAGRANGE", "JACOBI_20_00" "LAGRANGE", "JACOBI_20_00" "LAGRANGE", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" "CARTESIAN" "CARTESIAN" 
    Approximation Orders="FIRST", "THIRD" "FIRST", "THIRD" "FIRST", "THIRD" 
    n_dofs()=5535
    n_local_dofs()=1044
    n_constrained_dofs()=135
    n_local_constrained_dofs()=135
    n_vectors()=1
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 59.5345
      Average Off-Processor Bandwidth <= 2.96552
      Maximum  On-Processor Bandwidth <= 81
      Maximum Off-Processor Bandwidth <= 36
    DofMap Constraints
      Number of DoF Constraints = 135
      Average DoF Constraint Length= 0
      Number of Node Constraints = 0
   System #1, "StressSystem"
    Type "Explicit"
    Variables="sigma_00" "sigma_01" "sigma_02" "sigma_10" "sigma_11" "sigma_12" "sigma_20" "sigma_21" "sigma_22" "vonMises" 
    Finite Element Types="MONOMIAL", "JACOBI_20_00" "MONOMIAL", "JACOBI_20_00" "MONOMIAL", "JACOBI_20_00" "MONOMIAL", "JACOBI_20_00" "MONOMIAL", "JACOBI_20_00" "MONOMIAL", "JACOBI_20_00" "MONOMIAL", "JACOBI_20_00" "MONOMIAL", "JACOBI_20_00" "MONOMIAL", "JACOBI_20_00" "MONOMIAL", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" "CARTESIAN" "CARTESIAN" "CARTESIAN" "CARTESIAN" "CARTESIAN" "CARTESIAN" "CARTESIAN" "CARTESIAN" "CARTESIAN" 
    Approximation Orders="CONSTANT", "THIRD" "CONSTANT", "THIRD" "CONSTANT", "THIRD" "CONSTANT", "THIRD" "CONSTANT", "THIRD" "CONSTANT", "THIRD" "CONSTANT", "THIRD" "CONSTANT", "THIRD" "CONSTANT", "THIRD" "CONSTANT", "THIRD" 
    n_dofs()=12800
    n_local_dofs()=2140
    n_constrained_dofs()=0
    n_local_constrained_dofs()=0
    n_vectors()=1
    n_matrices()=0
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 0
      Average Off-Processor Bandwidth <= 0
      Maximum  On-Processor Bandwidth <= 0
      Maximum Off-Processor Bandwidth <= 0
    DofMap Constraints
      Number of DoF Constraints = 0
      Number of Node Constraints = 0

Initializing training parameters with random training set...
Parameter load_Fx: log scaling = 0
Parameter load_Fy: log scaling = 0
Parameter load_Fz: log scaling = 0
Parameter x_scaling: log scaling = 0


RBConstruction parameters:
system name: RBElasticity
constrained_problem: 0
Nmax: 15
Basis training error tolerance: 0.001
Aq operators attached: 3
Fq functions attached: 3
n_outputs: 0
Number of parameters: 4
Parameter load_Fx: Min = -5, Max = 5, value = 1
Parameter load_Fy: Min = -5, Max = 5, value = 1
Parameter load_Fz: Min = -5, Max = 5, value = 1
Parameter x_scaling: Min = 0.5, Max = 2, value = 1
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
load_Fx: -0.114724
load_Fy: 0.253877
load_Fz: -1.06179
x_scaling: 1.58763

Warning: Linear solver may not have converged! Final linear residual = 0.00367043, number of iterations = 5000

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 1 ----
Performing RB solves on training set
Maximum (relative) error bound is 132.46

Performing truth solve at parameter:
load_Fx: -3.82457
load_Fy: -1.37368
load_Fz: -0.0869978
x_scaling: 0.597728

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 2 ----
Performing RB solves on training set
Maximum (relative) error bound is 1.54107

Performing truth solve at parameter:
load_Fx: -3.48446
load_Fy: 0.70637
load_Fz: -0.0996004
x_scaling: 0.760771

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 3 ----
Performing RB solves on training set
Maximum (relative) error bound is 1.2238

Performing truth solve at parameter:
load_Fx: -1.68494
load_Fy: 0.33588
load_Fz: -0.0334245
x_scaling: 1.98747

Warning: Linear solver may not have converged! Final linear residual = 3.5854e-05, number of iterations = 5000

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 4 ----
Performing RB solves on training set
Maximum (relative) error bound is 0.293716

Performing truth solve at parameter:
load_Fx: -4.3869
load_Fy: 2.27439
load_Fz: 4.52502
x_scaling: 0.525225

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 5 ----
Performing RB solves on training set
Maximum (relative) error bound is 0.190118

Performing truth solve at parameter:
load_Fx: 4.98569
load_Fy: -0.199169
load_Fz: -0.333577
x_scaling: 0.737557

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 6 ----
Performing RB solves on training set
Maximum (relative) error bound is 0.0666646

Performing truth solve at parameter:
load_Fx: -4.13276
load_Fy: -2.96578
load_Fz: -2.05843
x_scaling: 1.90802

Warning: Linear solver may not have converged! Final linear residual = 3.07752e-05, number of iterations = 5000

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 7 ----
Performing RB solves on training set
Maximum (relative) error bound is 0.0370728

Performing truth solve at parameter:
load_Fx: 4.60009
load_Fy: 1.33364
load_Fz: 0.0860635
x_scaling: 1.30223

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 8 ----
Performing RB solves on training set
Maximum (relative) error bound is 0.0254904

Performing truth solve at parameter:
load_Fx: -3.22867
load_Fy: 0.860594
load_Fz: -0.0980799
x_scaling: 0.527005

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 9 ----
Performing RB solves on training set
Maximum (relative) error bound is 0.00165634

Performing truth solve at parameter:
load_Fx: -1.41026
load_Fy: -4.0215
load_Fz: 0.146317
x_scaling: 1.95283

Warning: Linear solver may not have converged! Final linear residual = 0.0111453, number of iterations = 5000

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 10 ----
Performing RB solves on training set
Maximum (relative) error bound is 0.00165601

Exiting greedy because the same parameters were selected twice
************************************************************************************************************************
***             WIDEN YOUR WINDOW TO 120 CHARACTERS.  Use 'enscript -r -fCourier9' to print this document            ***
************************************************************************************************************************

---------------------------------------------- PETSc Performance Summary: ----------------------------------------------

./reduced_basis_ex5-opt on a intel-11. named daedalus with 6 processors, by roystgnr Fri Aug 24 15:25:56 2012
Using Petsc Release Version 3.1.0, Patch 5, Mon Sep 27 11:51:54 CDT 2010

                         Max       Max/Min        Avg      Total 
Time (sec):           1.888e+02      1.00001   1.888e+02
Objects:              2.670e+02      1.00000   2.670e+02
Flops:                6.812e+10      1.78782   5.582e+10  3.349e+11
Flops/sec:            3.609e+08      1.78782   2.957e+08  1.774e+09
MPI Messages:         8.563e+04      2.00000   7.136e+04  4.282e+05
MPI Message Lengths:  1.076e+08      2.03925   1.213e+03  5.193e+08
MPI Reductions:       8.785e+04      1.00000

Flop counting convention: 1 flop = 1 real number operation of type (multiply/divide/add/subtract)
                            e.g., VecAXPY() for real vectors of length N --> 2N flops
                            and VecAXPY() for complex vectors of length N --> 8N flops

Summary of Stages:   ----- Time ------  ----- Flops -----  --- Messages ---  -- Message Lengths --  -- Reductions --
                        Avg     %Total     Avg     %Total   counts   %Total     Avg         %Total   counts   %Total 
 0:      Main Stage: 1.8877e+02 100.0%  3.3489e+11 100.0%  4.282e+05 100.0%  1.213e+03      100.0%  8.768e+04  99.8% 

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

VecDot              1091 1.0 1.6982e-01 1.9 2.28e+06 1.4 0.0e+00 0.0e+00 1.1e+03  0  0  0  0  1   0  0  0  0  1    71
VecMDot            40358 1.0 6.4002e+01 2.1 1.30e+09 1.4 0.0e+00 0.0e+00 4.0e+04 24  2  0  0 46  24  2  0  0 46   107
VecNorm            41766 1.0 5.0329e+01 2.8 8.72e+07 1.4 0.0e+00 0.0e+00 4.2e+04 16  0  0  0 48  16  0  0  0 48     9
VecScale           41793 1.0 8.6627e-02 1.6 4.36e+07 1.4 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0  2670
VecCopy             4256 1.0 9.9638e-03 1.4 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecSet             44592 1.0 9.7214e-02 1.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecAXPY             2808 1.0 1.7130e-02 2.8 5.86e+06 1.4 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0  1815
VecMAXPY           41723 1.0 7.0981e-01 1.5 1.38e+09 1.4 0.0e+00 0.0e+00 0.0e+00  0  2  0  0  0   0  2  0  0  0 10322
VecAssemblyBegin     270 1.0 2.6900e-01 2.1 0.00e+00 0.0 3.0e+01 2.6e+03 8.1e+02  0  0  0  0  1   0  0  0  0  1     0
VecAssemblyEnd       270 1.0 2.6155e-04 1.3 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecScatterBegin    42794 1.0 4.9415e-01 1.8 0.00e+00 0.0 4.3e+05 1.2e+03 0.0e+00  0  0100 99  0   0  0100 99  0     0
VecScatterEnd      42794 1.0 5.0049e+0110.3 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00 13  0  0  0  0  13  0  0  0  0     0
VecNormalize       41723 1.0 5.0407e+01 2.8 1.31e+08 1.4 0.0e+00 0.0e+00 4.2e+04 16  0  0  0 47  16  0  0  0 48    14
MatMult            41723 1.0 5.6691e+01 5.2 5.61e+09 1.4 4.2e+05 1.2e+03 0.0e+00 17  9 97 97  0  17  9 97 97  0   517
MatMultAdd          1028 1.0 1.5491e-01 2.0 1.39e+08 1.4 1.0e+04 1.2e+03 0.0e+00  0  0  2  2  0   0  0  2  2  0  4697
MatSolve           41723 1.0 8.2990e+01 1.8 5.66e+10 1.8 0.0e+00 0.0e+00 0.0e+00 35 82  0  0  0  35 82  0  0  0  3322
MatLUFactorNum        21 1.0 2.7763e+00 2.4 3.27e+09 2.4 0.0e+00 0.0e+00 0.0e+00  1  4  0  0  0   1  4  0  0  0  5118
MatILUFactorSym       21 1.0 1.0944e+01 2.3 0.00e+00 0.0 0.0e+00 0.0e+00 2.1e+01  5  0  0  0  0   5  0  0  0  0     0
MatAssemblyBegin    1181 1.0 6.4045e-01 2.9 0.00e+00 0.0 6.0e+01 5.6e+04 2.4e+03  0  0  0  1  3   0  0  0  1  3     0
MatAssemblyEnd      1181 1.0 2.5571e-01 1.2 0.00e+00 0.0 1.0e+02 3.0e+02 1.2e+03  0  0  0  0  1   0  0  0  0  1     0
MatGetRow          42804 1.4 1.9663e-02 1.8 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatGetRowIJ           21 1.0 2.5749e-05 8.3 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatGetOrdering        21 1.0 3.2234e-04 1.2 0.00e+00 0.0 0.0e+00 0.0e+00 4.2e+01  0  0  0  0  0   0  0  0  0  0     0
MatZeroEntries        35 1.0 6.0217e-03 1.9 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPGMRESOrthog     40358 1.0 6.4747e+01 2.1 2.59e+09 1.4 0.0e+00 0.0e+00 4.0e+04 25  4  0  0 46  25  4  0  0 46   212
KSPSetup              64 1.0 1.0443e-04 1.6 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPSolve              43 1.0 1.8655e+02 1.0 6.80e+10 1.8 4.2e+05 1.2e+03 8.2e+04 99100 97 97 94  99100 97 97 94  1791
PCSetUp               42 1.0 1.3733e+01 2.3 3.27e+09 2.4 0.0e+00 0.0e+00 6.3e+01  6  4  0  0  0   6  4  0  0  0  1035
PCSetUpOnBlocks       43 1.0 1.3732e+01 2.3 3.27e+09 2.4 0.0e+00 0.0e+00 6.3e+01  6  4  0  0  0   6  4  0  0  0  1035
PCApply            41723 1.0 8.3660e+01 1.7 5.66e+10 1.8 0.0e+00 0.0e+00 0.0e+00 35 82  0  0  0  35 82  0  0  0  3296
------------------------------------------------------------------------------------------------------------------------

Memory usage is given in bytes:

Object Type          Creations   Destructions     Memory  Descendants' Mem.
Reports information only for process 0.

--- Event Stage 0: Main Stage

                 Vec   137            137      1202848     0
         Vec Scatter     8              8         6944     0
           Index Set    79             79       309416     0
   IS L to G Mapping     2              2        15688     0
              Matrix    37             37    175148716     0
       Krylov Solver     2              2        18880     0
      Preconditioner     2              2         1408     0
========================================================================================================================
Average time to get PetscTime(): 0
Average time for MPI_Barrier(): 4.01974e-05
Average time for zero size MPI_Send(): 4.63327e-05
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
| Time:           Fri Aug 24 15:25:56 2012                         |
| OS:             Linux                                            |
| HostName:       daedalus                                         |
| OS Release:     2.6.32-34-generic                                |
| OS Version:     #76-Ubuntu SMP Tue Aug 30 17:05:01 UTC 2011      |
| Machine:        x86_64                                           |
| Username:       roystgnr                                         |
| Configuration:  ./configure run on Wed Aug 22 12:44:06 CDT 2012  |
-------------------------------------------------------------------
 -------------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=188.988, Active time=188.72                                                       |
 -------------------------------------------------------------------------------------------------------------------
| Event                                 nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                                 w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|-------------------------------------------------------------------------------------------------------------------|
|                                                                                                                   |
|                                                                                                                   |
| DofMap                                                                                                            |
|   add_neighbors_to_send_list()        2         0.0004      0.000192    0.0006      0.000276    0.00     0.00     |
|   build_constraint_matrix()           1498      0.0037      0.000002    0.0037      0.000002    0.00     0.00     |
|   build_sparsity()                    1         0.0087      0.008670    0.0096      0.009637    0.00     0.01     |
|   cnstrn_elem_mat_vec()               1498      0.0110      0.000007    0.0110      0.000007    0.01     0.01     |
|   create_dof_constraints()            2         0.0056      0.002823    0.0072      0.003602    0.00     0.00     |
|   distribute_dofs()                   2         0.0027      0.001339    0.0258      0.012903    0.00     0.01     |
|   dof_indices()                       11266     0.0057      0.000001    0.0057      0.000001    0.00     0.00     |
|   prepare_send_list()                 2         0.0001      0.000028    0.0001      0.000028    0.00     0.00     |
|   reinit()                            2         0.0072      0.003576    0.0072      0.003576    0.00     0.00     |
|                                                                                                                   |
| FE                                                                                                                |
|   compute_shape_functions()           5642      0.0301      0.000005    0.0301      0.000005    0.02     0.02     |
|   init_shape_functions()              2660      0.0183      0.000007    0.0183      0.000007    0.01     0.01     |
|                                                                                                                   |
| FEMap                                                                                                             |
|   compute_affine_map()                5642      0.0070      0.000001    0.0070      0.000001    0.00     0.00     |
|   compute_face_map()                  2646      0.0037      0.000001    0.0037      0.000001    0.00     0.00     |
|   init_face_shape_functions()         14        0.0001      0.000005    0.0001      0.000005    0.00     0.00     |
|   init_reference_to_physical_map()    2660      0.0134      0.000005    0.0134      0.000005    0.01     0.01     |
|                                                                                                                   |
| Mesh                                                                                                              |
|   find_neighbors()                    1         0.0020      0.002048    0.0023      0.002326    0.00     0.00     |
|   renumber_nodes_and_elem()           2         0.0002      0.000081    0.0002      0.000081    0.00     0.00     |
|                                                                                                                   |
| MeshCommunication                                                                                                 |
|   assign_global_indices()             1         0.0124      0.012448    0.0245      0.024546    0.01     0.01     |
|   compute_hilbert_indices()           2         0.0043      0.002133    0.0043      0.002133    0.00     0.00     |
|   find_global_indices()               2         0.0026      0.001322    0.0132      0.006611    0.00     0.01     |
|   parallel_sort()                     2         0.0013      0.000632    0.0053      0.002652    0.00     0.00     |
|                                                                                                                   |
| MeshTools::Generation                                                                                             |
|   build_cube()                        1         0.0007      0.000657    0.0007      0.000657    0.00     0.00     |
|                                                                                                                   |
| MetisPartitioner                                                                                                  |
|   partition()                         1         0.0028      0.002762    0.0103      0.010321    0.00     0.01     |
|                                                                                                                   |
| Parallel                                                                                                          |
|   allgather()                         16        0.0160      0.000999    0.0160      0.000999    0.01     0.01     |
|   barrier()                           11        0.0001      0.000005    0.0001      0.000005    0.00     0.00     |
|   broadcast()                         23        0.0012      0.000050    0.0012      0.000050    0.00     0.00     |
|   gather()                            121       0.0037      0.000030    0.0037      0.000030    0.00     0.00     |
|   max(scalar)                         3         0.0023      0.000776    0.0023      0.000776    0.00     0.00     |
|   max(vector)                         3         0.0001      0.000044    0.0001      0.000044    0.00     0.00     |
|   maxloc(scalar)                      11        0.0837      0.007605    0.0837      0.007605    0.04     0.04     |
|   min(vector)                         3         0.0009      0.000302    0.0009      0.000302    0.00     0.00     |
|   probe()                             90        0.0121      0.000135    0.0121      0.000135    0.01     0.01     |
|   receive()                           810       0.0005      0.000001    0.0126      0.000016    0.00     0.01     |
|   send()                              210       0.0003      0.000001    0.0003      0.000001    0.00     0.00     |
|   send_receive()                      98        0.0002      0.000002    0.0128      0.000130    0.00     0.01     |
|   sum()                               25        0.0074      0.000297    0.0074      0.000297    0.00     0.00     |
|                                                                                                                   |
| Parallel::Request                                                                                                 |
|   wait()                              810       0.0004      0.000001    0.0004      0.000001    0.00     0.00     |
|                                                                                                                   |
| Partitioner                                                                                                       |
|   set_node_processor_ids()            1         0.0002      0.000208    0.0011      0.001121    0.00     0.00     |
|   set_parent_processor_ids()          1         0.0001      0.000078    0.0001      0.000078    0.00     0.00     |
|                                                                                                                   |
| PetscLinearSolver                                                                                                 |
|   solve()                             43        186.6042    4.339634    186.6042    4.339634    98.88    98.88    |
|                                                                                                                   |
| RBConstruction                                                                                                    |
|   add_scaled_matrix_and_vector()      7         0.2908      0.041540    0.3835      0.054785    0.15     0.20     |
|   clear()                             1         0.0008      0.000837    0.0008      0.000837    0.00     0.00     |
|   compute_Fq_representor_innerprods() 1         0.5963      0.596333    2.3541      2.354120    0.32     1.25     |
|   compute_max_error_bound()           11        0.0113      0.001025    0.2251      0.020464    0.01     0.12     |
|   enrich_RB_space()                   10        0.0304      0.003042    0.0304      0.003042    0.02     0.02     |
|   train_reduced_basis()               1         0.0199      0.019899    188.2090    188.209041  0.01     99.73    |
|   truth_assembly()                    10        0.1666      0.016658    0.1666      0.016658    0.09     0.09     |
|   truth_solve()                       10        0.0209      0.002092    169.8593    16.985930   0.01     90.01    |
|   update_RB_system_matrices()         10        0.1469      0.014686    0.1469      0.014686    0.08     0.08     |
|   update_residual_terms()             10        0.3986      0.039861    15.5733     1.557328    0.21     8.25     |
|                                                                                                                   |
| RBEvaluation                                                                                                      |
|   clear()                             1         0.0001      0.000064    0.0001      0.000064    0.00     0.00     |
|   compute_residual_dual_norm()        1837      0.1147      0.000062    0.1147      0.000062    0.06     0.06     |
|   rb_solve()                          1837      0.0120      0.000007    0.1269      0.000069    0.01     0.07     |
|   resize_data_structures()            1         0.0000      0.000045    0.0000      0.000045    0.00     0.00     |
|   write_offline_data_to_files()       1         0.0339      0.033867    0.0339      0.033867    0.02     0.02     |
 -------------------------------------------------------------------------------------------------------------------
| Totals:                               39576     188.7203                                        100.00            |
 -------------------------------------------------------------------------------------------------------------------

 Mesh Information:
  mesh_dimension()=3
  spatial_dimension()=3
  n_nodes()=1845
    n_local_nodes()=348
  n_elem()=1280
    n_local_elem()=214
    n_active_elem()=1280
  n_subdomains()=1
  n_partitions()=6
  n_processors()=6
  n_threads()=1
  processor_id()=0

*** Warning, This code is untested, experimental, or likely to see future API changes: src/reduced_basis/rb_parametrized.C, line 40, compiled Aug 24 2012 at 15:15:42 ***
 EquationSystems
  n_systems()=2
   System #0, "RBElasticity"
    Type "RBConstruction"
    Variables="u" "v" "w" 
    Finite Element Types="LAGRANGE", "JACOBI_20_00" "LAGRANGE", "JACOBI_20_00" "LAGRANGE", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" "CARTESIAN" "CARTESIAN" 
    Approximation Orders="FIRST", "THIRD" "FIRST", "THIRD" "FIRST", "THIRD" 
    n_dofs()=5535
    n_local_dofs()=1044
    n_constrained_dofs()=135
    n_local_constrained_dofs()=135
    n_vectors()=1
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 59.5345
      Average Off-Processor Bandwidth <= 2.96552
      Maximum  On-Processor Bandwidth <= 81
      Maximum Off-Processor Bandwidth <= 36
    DofMap Constraints
      Number of DoF Constraints = 135
      Average DoF Constraint Length= 0
      Number of Node Constraints = 0
   System #1, "StressSystem"
    Type "Explicit"
    Variables="sigma_00" "sigma_01" "sigma_02" "sigma_10" "sigma_11" "sigma_12" "sigma_20" "sigma_21" "sigma_22" "vonMises" 
    Finite Element Types="MONOMIAL", "JACOBI_20_00" "MONOMIAL", "JACOBI_20_00" "MONOMIAL", "JACOBI_20_00" "MONOMIAL", "JACOBI_20_00" "MONOMIAL", "JACOBI_20_00" "MONOMIAL", "JACOBI_20_00" "MONOMIAL", "JACOBI_20_00" "MONOMIAL", "JACOBI_20_00" "MONOMIAL", "JACOBI_20_00" "MONOMIAL", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" "CARTESIAN" "CARTESIAN" "CARTESIAN" "CARTESIAN" "CARTESIAN" "CARTESIAN" "CARTESIAN" "CARTESIAN" "CARTESIAN" 
    Approximation Orders="CONSTANT", "THIRD" "CONSTANT", "THIRD" "CONSTANT", "THIRD" "CONSTANT", "THIRD" "CONSTANT", "THIRD" "CONSTANT", "THIRD" "CONSTANT", "THIRD" "CONSTANT", "THIRD" "CONSTANT", "THIRD" "CONSTANT", "THIRD" 
    n_dofs()=12800
    n_local_dofs()=2140
    n_constrained_dofs()=0
    n_local_constrained_dofs()=0
    n_vectors()=1
    n_matrices()=0
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 0
      Average Off-Processor Bandwidth <= 0
      Maximum  On-Processor Bandwidth <= 0
      Maximum Off-Processor Bandwidth <= 0
    DofMap Constraints
      Number of DoF Constraints = 0
      Number of Node Constraints = 0

load_Fx: 0
load_Fy: 0
load_Fz: -1
x_scaling: 1.3

************************************************************************************************************************
***             WIDEN YOUR WINDOW TO 120 CHARACTERS.  Use 'enscript -r -fCourier9' to print this document            ***
************************************************************************************************************************

---------------------------------------------- PETSc Performance Summary: ----------------------------------------------

./reduced_basis_ex5-opt on a intel-11. named daedalus with 6 processors, by roystgnr Fri Aug 24 15:25:57 2012
Using Petsc Release Version 3.1.0, Patch 5, Mon Sep 27 11:51:54 CDT 2010

                         Max       Max/Min        Avg      Total 
Time (sec):           2.197e-01      1.11407   2.014e-01
Objects:              3.000e+01      1.00000   3.000e+01
Flops:                2.088e+04      1.39200   1.845e+04  1.107e+05
Flops/sec:            1.052e+05      1.38308   9.156e+04  5.494e+05
MPI Messages:         1.200e+01      2.00000   1.000e+01  6.000e+01
MPI Message Lengths:  1.481e+04      2.05582   1.183e+03  7.098e+04
MPI Reductions:       7.800e+01      1.00000

Flop counting convention: 1 flop = 1 real number operation of type (multiply/divide/add/subtract)
                            e.g., VecAXPY() for real vectors of length N --> 2N flops
                            and VecAXPY() for complex vectors of length N --> 8N flops

Summary of Stages:   ----- Time ------  ----- Flops -----  --- Messages ---  -- Message Lengths --  -- Reductions --
                        Avg     %Total     Avg     %Total   counts   %Total     Avg         %Total   counts   %Total 
 0:      Main Stage: 2.0130e-01 100.0%  1.1070e+05 100.0%  6.000e+01 100.0%  1.183e+03      100.0%  5.300e+01  67.9% 

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

VecCopy                2 1.0 1.7881e-05 3.4 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecSet                18 1.0 8.0585e-05 8.7 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecAXPY               10 1.0 6.5804e-05 2.0 2.09e+04 1.4 0.0e+00 0.0e+00 0.0e+00  0100  0  0  0   0100  0  0  0  1682
VecAssemblyBegin      13 1.0 1.4096e-02 6.9 0.00e+00 0.0 0.0e+00 0.0e+00 3.9e+01  4  0  0  0 50   4  0  0  0 74     0
VecAssemblyEnd        13 1.0 3.3617e-05 1.6 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecScatterBegin        2 1.0 8.3113e-0426.8 0.00e+00 0.0 2.0e+01 2.4e+03 0.0e+00  0  0 33 67  0   0  0 33 67  0     0
VecScatterEnd          2 1.0 8.3280e-0491.9 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatZeroEntries         2 1.0 5.3811e-04 2.4 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
------------------------------------------------------------------------------------------------------------------------

Memory usage is given in bytes:

Object Type          Creations   Destructions     Memory  Descendants' Mem.
Reports information only for process 0.

--- Event Stage 0: Main Stage

                 Vec    19             19       197984     0
         Vec Scatter     2              2         1736     0
           Index Set     4              4         4208     0
   IS L to G Mapping     2              2        15688     0
              Matrix     3              3       814592     0
========================================================================================================================
Average time to get PetscTime(): 9.53674e-08
Average time for MPI_Barrier(): 3.60012e-05
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
| Time:           Fri Aug 24 15:25:57 2012                         |
| OS:             Linux                                            |
| HostName:       daedalus                                         |
| OS Release:     2.6.32-34-generic                                |
| OS Version:     #76-Ubuntu SMP Tue Aug 30 17:05:01 UTC 2011      |
| Machine:        x86_64                                           |
| Username:       roystgnr                                         |
| Configuration:  ./configure run on Wed Aug 22 12:44:06 CDT 2012  |
-------------------------------------------------------------------
 ----------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=0.394825, Active time=0.171761                                                 |
 ----------------------------------------------------------------------------------------------------------------
| Event                              nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                              w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|----------------------------------------------------------------------------------------------------------------|
|                                                                                                                |
|                                                                                                                |
| DofMap                                                                                                         |
|   add_neighbors_to_send_list()     2         0.0004      0.000212    0.0006      0.000317    0.25     0.37     |
|   build_sparsity()                 1         0.0093      0.009295    0.0106      0.010594    5.41     6.17     |
|   create_dof_constraints()         2         0.0074      0.003713    0.0095      0.004745    4.32     5.53     |
|   distribute_dofs()                2         0.0029      0.001460    0.0288      0.014378    1.70     16.74    |
|   dof_indices()                    10838     0.0050      0.000000    0.0050      0.000000    2.92     2.92     |
|   prepare_send_list()              2         0.0001      0.000033    0.0001      0.000033    0.04     0.04     |
|   reinit()                         2         0.0076      0.003818    0.0076      0.003818    4.45     4.45     |
|                                                                                                                |
| EquationSystems                                                                                                |
|   build_solution_vector()          1         0.0027      0.002676    0.0112      0.011200    1.56     6.52     |
|                                                                                                                |
| ExodusII_IO                                                                                                    |
|   write_nodal_data()               1         0.0029      0.002893    0.0029      0.002893    1.68     1.68     |
|                                                                                                                |
| FE                                                                                                             |
|   compute_shape_functions()        214       0.0003      0.000001    0.0003      0.000001    0.17     0.17     |
|   init_shape_functions()           1         0.0000      0.000020    0.0000      0.000020    0.01     0.01     |
|                                                                                                                |
| FEMap                                                                                                          |
|   compute_affine_map()             214       0.0003      0.000001    0.0003      0.000001    0.15     0.15     |
|   init_reference_to_physical_map() 1         0.0000      0.000031    0.0000      0.000031    0.02     0.02     |
|                                                                                                                |
| Mesh                                                                                                           |
|   find_neighbors()                 1         0.0020      0.001962    0.0058      0.005822    1.14     3.39     |
|   renumber_nodes_and_elem()        2         0.0002      0.000095    0.0002      0.000095    0.11     0.11     |
|                                                                                                                |
| MeshCommunication                                                                                              |
|   assign_global_indices()          1         0.0113      0.011315    0.0322      0.032172    6.59     18.73    |
|   compute_hilbert_indices()        2         0.0070      0.003490    0.0070      0.003490    4.06     4.06     |
|   find_global_indices()            2         0.0006      0.000281    0.0131      0.006552    0.33     7.63     |
|   parallel_sort()                  2         0.0045      0.002263    0.0046      0.002302    2.64     2.68     |
|                                                                                                                |
| MeshOutput                                                                                                     |
|   write_equation_systems()         1         0.0000      0.000031    0.0141      0.014124    0.02     8.22     |
|                                                                                                                |
| MeshTools::Generation                                                                                          |
|   build_cube()                     1         0.0007      0.000675    0.0007      0.000675    0.39     0.39     |
|                                                                                                                |
| MetisPartitioner                                                                                               |
|   partition()                      1         0.0033      0.003334    0.0107      0.010673    1.94     6.21     |
|                                                                                                                |
| Parallel                                                                                                       |
|   allgather()                      76        0.0422      0.000555    0.0422      0.000555    24.57    24.57    |
|   barrier()                        1         0.0000      0.000017    0.0000      0.000017    0.01     0.01     |
|   broadcast()                      466       0.0026      0.000006    0.0026      0.000006    1.54     1.53     |
|   gather()                         1         0.0003      0.000267    0.0003      0.000267    0.16     0.16     |
|   max(scalar)                      3         0.0053      0.001760    0.0053      0.001760    3.07     3.07     |
|   max(vector)                      3         0.0007      0.000234    0.0007      0.000234    0.41     0.41     |
|   min(vector)                      3         0.0009      0.000314    0.0009      0.000314    0.55     0.55     |
|   probe()                          90        0.0071      0.000079    0.0071      0.000079    4.15     4.15     |
|   receive()                        90        0.0003      0.000003    0.0074      0.000082    0.16     4.32     |
|   send()                           90        0.0002      0.000002    0.0002      0.000002    0.11     0.11     |
|   send_receive()                   98        0.0002      0.000002    0.0079      0.000081    0.13     4.61     |
|   sum()                            18        0.0097      0.000537    0.0097      0.000537    5.63     5.63     |
|                                                                                                                |
| Parallel::Request                                                                                              |
|   wait()                           90        0.0001      0.000001    0.0001      0.000001    0.04     0.04     |
|                                                                                                                |
| Partitioner                                                                                                    |
|   set_node_processor_ids()         1         0.0003      0.000264    0.0015      0.001508    0.15     0.88     |
|   set_parent_processor_ids()       1         0.0001      0.000098    0.0001      0.000098    0.06     0.06     |
|                                                                                                                |
| RBConstruction                                                                                                 |
|   clear()                          1         0.0002      0.000153    0.0002      0.000153    0.09     0.09     |
|   load_rb_solution()               1         0.0003      0.000329    0.0003      0.000329    0.19     0.19     |
|                                                                                                                |
| RBEvaluation                                                                                                   |
|   clear()                          1         0.0001      0.000055    0.0001      0.000055    0.03     0.03     |
|   compute_residual_dual_norm()     1         0.0001      0.000148    0.0001      0.000148    0.09     0.09     |
|   rb_solve()                       1         0.0111      0.011140    0.0113      0.011289    6.49     6.57     |
|   read_offline_data_from_files()   1         0.0008      0.000799    0.0008      0.000827    0.47     0.48     |
|   resize_data_structures()         1         0.0000      0.000027    0.0000      0.000027    0.02     0.02     |
|                                                                                                                |
| main                                                                                                           |
|   compute_stresses()               1         0.0206      0.020617    0.0223      0.022310    12.00    12.99    |
 ----------------------------------------------------------------------------------------------------------------
| Totals:                            12333     0.1718                                          100.00            |
 ----------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running  mpirun -np 6 ./reduced_basis_ex5-opt -ksp_type cg -online_mode ?
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
