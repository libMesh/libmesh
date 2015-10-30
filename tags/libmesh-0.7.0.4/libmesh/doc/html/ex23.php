<?php $root=""; ?>
<?php require($root."navigation.php"); ?>
<html>
<head>
  <?php load_style($root); ?>
</head>
 
<body>
 
<?php make_navigation("ex23",$root)?>
 
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
        #include "gmv_io.h"
        #include "equation_systems.h"
        #include "fe.h"
        #include "quadrature_gauss.h"
        #include "dof_map.h"
        #include "sparse_matrix.h"
        #include "numeric_vector.h"
        #include "dense_matrix.h"
        #include "dense_vector.h"
        #include "fe_interface.h"
        #include "getpot.h"
        #include "o_string_stream.h"
        #include "elem.h"
        #include "simple_rb.h"
        
</pre>
</div>
<div class = "comment">
In this example problem we use the Certified Reduced Basis method
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
specified in ex23.in.


<br><br>We also attach four outputs to the system which are averages over certain
subregions of the domain. In Online mode, we print out the values of these
outputs as well as rigorous error bounds with respect to the output
associated with the "truth" finite element discretization.


<br><br>

<br><br>Populate the sets non-Dirichlet and Dirichlet dofs,
used to impose boundary conditions
</div>

<div class ="fragment">
<pre>
        void initialize_dirichlet_dofs(FEMContext &c,
                                       System& system,
                                       std::set&lt;unsigned int&gt;& dirichlet_dofs_set);
        
</pre>
</div>
<div class = "comment">
Expansion of the PDE operator
</div>

<div class ="fragment">
<pre>
        Number theta_a_0(RBThetaData& ) { return 0.05; }
        void A0(FEMContext&, System&);
        Number theta_a_1(RBThetaData& data) { return data.get_mu()[0]; }
        void A1(FEMContext&, System&);
        Number theta_a_2(RBThetaData& data) { return data.get_mu()[1]; }
        void A2(FEMContext&, System&);
        
</pre>
</div>
<div class = "comment">
Expansion of the RHS
</div>

<div class ="fragment">
<pre>
        Number theta_f_0(RBThetaData& ) { return 1.; }
        void F0(FEMContext&, System&);
        
</pre>
</div>
<div class = "comment">
Expansion of the outputs
</div>

<div class ="fragment">
<pre>
        Number theta_l(RBThetaData& ) { return 1.; }
        void L0_assembly(FEMContext&, System&);
        void L1_assembly(FEMContext&, System&);
        void L2_assembly(FEMContext&, System&);
        void L3_assembly(FEMContext&, System&);
        
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
        #elif !defined(LIBMESH_HAVE_PETSC)
</pre>
</div>
<div class = "comment">
FIXME: This example currently segfaults with Trilinos?
</div>

<div class ="fragment">
<pre>
          libmesh_example_assert(false, "--enable-petsc");
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
Parse the input file (ex23.in) using GetPot
</div>

<div class ="fragment">
<pre>
          std::string parameters_filename = "ex23.in";
          GetPot infile(parameters_filename);
        
          unsigned int n_elem = infile("n_elem", 1);       // Determines the number of elements in the "truth" mesh
          const unsigned int dim = 2;                      // The number of spatial dimensions
          
          bool online_mode = infile("online_mode", false); // Are we in Online mode?
        
</pre>
</div>
<div class = "comment">
Create a one-dimensional mesh.
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
We override RBSystem with SimpleRB in order to specialize a few functions.
</div>

<div class ="fragment">
<pre>
          SimpleRB & system =
            equation_systems.add_system&lt;SimpleRB&gt; ("RBConvectionDiffusion");
        
</pre>
</div>
<div class = "comment">
Point the systems to the input file defining the problem
</div>

<div class ="fragment">
<pre>
          system.parameters_filename = parameters_filename;
        
          system.attach_dirichlet_dof_initialization(initialize_dirichlet_dofs);
          
</pre>
</div>
<div class = "comment">
Attach the expansion of the PDE operator. The third argument
here refers to assembly over the boundary of the mesh, but this
problem only requires internal assembly and hence it is set to NULL
</div>

<div class ="fragment">
<pre>
          system.attach_A_q(theta_a_0, A0, NULL);
          system.attach_A_q(theta_a_1, A1, NULL);
          system.attach_A_q(theta_a_2, A2, NULL);
          
</pre>
</div>
<div class = "comment">
Attach the expansion of the RHS
</div>

<div class ="fragment">
<pre>
          system.attach_F_q(theta_f_0, F0, NULL);
          
</pre>
</div>
<div class = "comment">
Attach output assembly
</div>

<div class ="fragment">
<pre>
          system.attach_output(theta_l, L0_assembly, NULL);
          system.attach_output(theta_l, L1_assembly, NULL);
          system.attach_output(theta_l, L2_assembly, NULL);
          system.attach_output(theta_l, L3_assembly, NULL);
          
</pre>
</div>
<div class = "comment">
We reuse the operator A0 as the inner product matrix
</div>

<div class ="fragment">
<pre>
          system.attach_inner_prod_assembly(A0);
        
</pre>
</div>
<div class = "comment">
Initialize the data structures for the equation system.
</div>

<div class ="fragment">
<pre>
          equation_systems.init ();
        
          if(system.initialize_calN_dependent_data)
          {
</pre>
</div>
<div class = "comment">
Print out some information about the "truth" discretization
</div>

<div class ="fragment">
<pre>
            mesh.print_info();
            equation_systems.print_info();
          }
          
</pre>
</div>
<div class = "comment">
Initialize the RB data structures.
If we're in Offline Mode (online_mode == false) then
also pre-assemble the Dirichlet dofs list, RB operators etc
</div>

<div class ="fragment">
<pre>
          system.initialize_RB_system(online_mode);
        
          if(!online_mode)
          {
</pre>
</div>
<div class = "comment">
Compute the reduced basis space by computing "snapshots", i.e.
"truth" solves, at well-chosen parameter values and employing
these snapshots as basis functions.
</div>

<div class ="fragment">
<pre>
            system.train_reduced_basis();
            
</pre>
</div>
<div class = "comment">
Write out the reduced basis data
</div>

<div class ="fragment">
<pre>
            system.write_offline_data_to_files();
          }
          else
          {
</pre>
</div>
<div class = "comment">
Read in the reduced basis data
</div>

<div class ="fragment">
<pre>
            system.read_offline_data_from_files();
            
</pre>
</div>
<div class = "comment">
Get the parameters at which we do a reduced basis solve
</div>

<div class ="fragment">
<pre>
            unsigned int online_N = infile("online_N",1);
            std::vector&lt;Real&gt; online_mu_vector(system.get_n_params());
            for(unsigned int i=0; i&lt;system.get_n_params(); i++)
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
            system.set_current_parameters(online_mu_vector);
            system.print_current_parameters();
        
</pre>
</div>
<div class = "comment">
Now do the Online solve using the precomputed reduced basis
</div>

<div class ="fragment">
<pre>
            system.RB_solve(online_N);
        
</pre>
</div>
<div class = "comment">
Print out outputs as well as the corresponding output error bounds.
</div>

<div class ="fragment">
<pre>
            std::cout &lt;&lt; "output 1, value = " &lt;&lt; system.RB_outputs[0]
                      &lt;&lt; ", bound = " &lt;&lt; system.RB_output_error_bounds[0]
                      &lt;&lt; std::endl;
            std::cout &lt;&lt; "output 2, value = " &lt;&lt; system.RB_outputs[1]
                      &lt;&lt; ", bound = " &lt;&lt; system.RB_output_error_bounds[1]
                      &lt;&lt; std::endl;
            std::cout &lt;&lt; "output 3, value = " &lt;&lt; system.RB_outputs[2]
                      &lt;&lt; ", bound = " &lt;&lt; system.RB_output_error_bounds[2]
                      &lt;&lt; std::endl;
            std::cout &lt;&lt; "output 4, value = " &lt;&lt; system.RB_outputs[3]
                      &lt;&lt; ", bound = " &lt;&lt; system.RB_output_error_bounds[3]
                      &lt;&lt; std::endl;
        
</pre>
</div>
<div class = "comment">
If we stored the basis functions in the offline_data directory, plot the RB solution
</div>

<div class ="fragment">
<pre>
            if(system.store_basis_functions)
            {
              system.load_RB_solution();
              GMVIO(mesh).write_equation_systems ("RB_sol.gmv",equation_systems);
              
              system.load_basis_function(0);
              GMVIO(mesh).write_equation_systems ("bf0.gmv",equation_systems);
            }
          }
        
          return 0;
        }
        
</pre>
</div>
<div class = "comment">
The Laplacian
</div>

<div class ="fragment">
<pre>
        void A0(FEMContext &c, System& system)
        {
          SimpleRB& rb_system = libmesh_cast_ref&lt;SimpleRB&&gt;(system);
          const unsigned int u_var = rb_system.u_var;
        
          const std::vector&lt;Real&gt; &JxW =
            c.element_fe_var[u_var]-&gt;get_JxW();
        
</pre>
</div>
<div class = "comment">
The velocity shape function gradients at interior
quadrature points.
</div>

<div class ="fragment">
<pre>
          const std::vector&lt;std::vector&lt;RealGradient&gt; &gt;& dphi =
            c.element_fe_var[u_var]-&gt;get_dphi();
        
</pre>
</div>
<div class = "comment">
The number of local degrees of freedom in each variable
</div>

<div class ="fragment">
<pre>
          const unsigned int n_u_dofs = c.dof_indices_var[u_var].size();
        
</pre>
</div>
<div class = "comment">
Now we will build the affine operator
</div>

<div class ="fragment">
<pre>
          unsigned int n_qpoints = c.element_qrule-&gt;n_points();
        
          for (unsigned int qp=0; qp != n_qpoints; qp++)
            for (unsigned int i=0; i != n_u_dofs; i++)
              for (unsigned int j=0; j != n_u_dofs; j++)
                c.elem_jacobian(i,j) += JxW[qp] * dphi[j][qp]*dphi[i][qp];
        }
        
</pre>
</div>
<div class = "comment">
Convection in the x-direction
</div>

<div class ="fragment">
<pre>
        void A1(FEMContext &c, System& system)
        {
          SimpleRB& rb_system = libmesh_cast_ref&lt;SimpleRB&&gt;(system);
          const unsigned int u_var = rb_system.u_var;
        
          const std::vector&lt;Real&gt; &JxW =
            c.element_fe_var[u_var]-&gt;get_JxW();
        
          const std::vector&lt;std::vector&lt;Real&gt; &gt;& phi =
            c.element_fe_var[u_var]-&gt;get_phi();
        
          const std::vector&lt;std::vector&lt;RealGradient&gt; &gt;& dphi =
            c.element_fe_var[u_var]-&gt;get_dphi();
        
</pre>
</div>
<div class = "comment">
The number of local degrees of freedom in each variable
</div>

<div class ="fragment">
<pre>
          const unsigned int n_u_dofs = c.dof_indices_var[u_var].size();
        
</pre>
</div>
<div class = "comment">
Now we will build the affine operator
</div>

<div class ="fragment">
<pre>
          unsigned int n_qpoints = c.element_qrule-&gt;n_points();
        
          for (unsigned int qp=0; qp != n_qpoints; qp++)
            for (unsigned int i=0; i != n_u_dofs; i++)
              for (unsigned int j=0; j != n_u_dofs; j++)
                c.elem_jacobian(i,j) += JxW[qp] *dphi[j][qp](0)*phi[i][qp];
        }
        
</pre>
</div>
<div class = "comment">
Convection in the y-direction
</div>

<div class ="fragment">
<pre>
        void A2(FEMContext &c, System& system)
        {
          SimpleRB& rb_system = libmesh_cast_ref&lt;SimpleRB&&gt;(system);
          const unsigned int u_var = rb_system.u_var;
        
          const std::vector&lt;Real&gt; &JxW =
            c.element_fe_var[u_var]-&gt;get_JxW();
        
          const std::vector&lt;std::vector&lt;Real&gt; &gt;& phi =
            c.element_fe_var[u_var]-&gt;get_phi();
        
          const std::vector&lt;std::vector&lt;RealGradient&gt; &gt;& dphi =
            c.element_fe_var[u_var]-&gt;get_dphi();
        
</pre>
</div>
<div class = "comment">
The number of local degrees of freedom in each variable
</div>

<div class ="fragment">
<pre>
          const unsigned int n_u_dofs = c.dof_indices_var[u_var].size();
        
</pre>
</div>
<div class = "comment">
Now we will build the affine operator
</div>

<div class ="fragment">
<pre>
          unsigned int n_qpoints = c.element_qrule-&gt;n_points();
        
          for (unsigned int qp=0; qp != n_qpoints; qp++)
            for (unsigned int i=0; i != n_u_dofs; i++)
              for (unsigned int j=0; j != n_u_dofs; j++)
                c.elem_jacobian(i,j) += JxW[qp] *dphi[j][qp](1)*phi[i][qp];
        }
        
</pre>
</div>
<div class = "comment">
Source term, 1 throughout the domain
</div>

<div class ="fragment">
<pre>
        void F0(FEMContext &c, System& system)
        {
          SimpleRB& rb_system = libmesh_cast_ref&lt;SimpleRB&&gt;(system);
          const unsigned int u_var = rb_system.u_var;
        
          const std::vector&lt;Real&gt; &JxW =
            c.element_fe_var[u_var]-&gt;get_JxW();
        
          const std::vector&lt;std::vector&lt;Real&gt; &gt;& phi =
            c.element_fe_var[u_var]-&gt;get_phi();
        
</pre>
</div>
<div class = "comment">
The number of local degrees of freedom in each variable
</div>

<div class ="fragment">
<pre>
          const unsigned int n_u_dofs = c.dof_indices_var[u_var].size();
        
</pre>
</div>
<div class = "comment">
Now we will build the affine operator
</div>

<div class ="fragment">
<pre>
          unsigned int n_qpoints = c.element_qrule-&gt;n_points();
        
          for (unsigned int qp=0; qp != n_qpoints; qp++)
            for (unsigned int i=0; i != n_u_dofs; i++)
              c.elem_residual(i) += JxW[qp] * ( 1.*phi[i][qp] );
        }
        
</pre>
</div>
<div class = "comment">
Output: Average value over the region [0.7,0.8]x[0.7,0.8]
</div>

<div class ="fragment">
<pre>
        void L0_assembly(FEMContext &c, System& system)
        {
          SimpleRB& rb_system = libmesh_cast_ref&lt;SimpleRB&&gt;(system);
          const unsigned int u_var = rb_system.u_var;
        
          const std::vector&lt;Real&gt; &JxW =
            c.element_fe_var[u_var]-&gt;get_JxW();
        
          const std::vector&lt;std::vector&lt;Real&gt; &gt;& phi =
            c.element_fe_var[u_var]-&gt;get_phi();
        
</pre>
</div>
<div class = "comment">
The number of local degrees of freedom in each variable
</div>

<div class ="fragment">
<pre>
          const unsigned int n_u_dofs = c.dof_indices_var[u_var].size();
        
</pre>
</div>
<div class = "comment">
Now we will build the affine operator
</div>

<div class ="fragment">
<pre>
          unsigned int n_qpoints = c.element_qrule-&gt;n_points();
        
          Point centroid = c.elem-&gt;centroid();
          if( (0.7 &lt;= centroid(0)) && (centroid(0) &lt;= 0.8) &&
              (0.7 &lt;= centroid(1)) && (centroid(1) &lt;= 0.8) )
            for (unsigned int qp=0; qp != n_qpoints; qp++)
              for (unsigned int i=0; i != n_u_dofs; i++)
                c.elem_residual(i) += JxW[qp] * ( 1.*phi[i][qp] ) / 0.01;
        }
        
</pre>
</div>
<div class = "comment">
Output: Average value over the region [0.2,0.3]x[0.7,0.8]
</div>

<div class ="fragment">
<pre>
        void L1_assembly(FEMContext &c, System& system)
        {
          SimpleRB& rb_system = libmesh_cast_ref&lt;SimpleRB&&gt;(system);
          const unsigned int u_var = rb_system.u_var;
        
          const std::vector&lt;Real&gt; &JxW =
            c.element_fe_var[u_var]-&gt;get_JxW();
        
          const std::vector&lt;std::vector&lt;Real&gt; &gt;& phi =
            c.element_fe_var[u_var]-&gt;get_phi();
        
</pre>
</div>
<div class = "comment">
The number of local degrees of freedom in each variable
</div>

<div class ="fragment">
<pre>
          const unsigned int n_u_dofs = c.dof_indices_var[u_var].size();
        
</pre>
</div>
<div class = "comment">
Now we will build the affine operator
</div>

<div class ="fragment">
<pre>
          unsigned int n_qpoints = c.element_qrule-&gt;n_points();
        
          Point centroid = c.elem-&gt;centroid();
          if( (0.2 &lt;= centroid(0)) && (centroid(0) &lt;= 0.3) &&
              (0.7 &lt;= centroid(1)) && (centroid(1) &lt;= 0.8) )
            for (unsigned int qp=0; qp != n_qpoints; qp++)
              for (unsigned int i=0; i != n_u_dofs; i++)
                c.elem_residual(i) += JxW[qp] * ( 1.*phi[i][qp] ) / 0.01;
        }
        
</pre>
</div>
<div class = "comment">
Output: Average value over the region [0.2,0.3]x[0.2,0.3]
</div>

<div class ="fragment">
<pre>
        void L2_assembly(FEMContext &c, System& system)
        {
          SimpleRB& rb_system = libmesh_cast_ref&lt;SimpleRB&&gt;(system);
          const unsigned int u_var = rb_system.u_var;
        
          const std::vector&lt;Real&gt; &JxW =
            c.element_fe_var[u_var]-&gt;get_JxW();
        
          const std::vector&lt;std::vector&lt;Real&gt; &gt;& phi =
            c.element_fe_var[u_var]-&gt;get_phi();
        
</pre>
</div>
<div class = "comment">
The number of local degrees of freedom in each variable
</div>

<div class ="fragment">
<pre>
          const unsigned int n_u_dofs = c.dof_indices_var[u_var].size();
        
</pre>
</div>
<div class = "comment">
Now we will build the affine operator
</div>

<div class ="fragment">
<pre>
          unsigned int n_qpoints = c.element_qrule-&gt;n_points();
        
          Point centroid = c.elem-&gt;centroid();
          if( (0.2 &lt;= centroid(0)) && (centroid(0) &lt;= 0.3) &&
              (0.2 &lt;= centroid(1)) && (centroid(1) &lt;= 0.3) )
            for (unsigned int qp=0; qp != n_qpoints; qp++)
              for (unsigned int i=0; i != n_u_dofs; i++)
                c.elem_residual(i) += JxW[qp] * ( 1.*phi[i][qp] ) / 0.01;
        }
        
</pre>
</div>
<div class = "comment">
Output: Average value over the region [0.7,0.8]x[0.2,0.3]
</div>

<div class ="fragment">
<pre>
        void L3_assembly(FEMContext &c, System& system)
        {
          SimpleRB& rb_system = libmesh_cast_ref&lt;SimpleRB&&gt;(system);
          const unsigned int u_var = rb_system.u_var;
        
          const std::vector&lt;Real&gt; &JxW =
            c.element_fe_var[u_var]-&gt;get_JxW();
        
          const std::vector&lt;std::vector&lt;Real&gt; &gt;& phi =
            c.element_fe_var[u_var]-&gt;get_phi();
        
</pre>
</div>
<div class = "comment">
The number of local degrees of freedom in each variable
</div>

<div class ="fragment">
<pre>
          const unsigned int n_u_dofs = c.dof_indices_var[u_var].size();
        
</pre>
</div>
<div class = "comment">
Now we will build the affine operator
</div>

<div class ="fragment">
<pre>
          unsigned int n_qpoints = c.element_qrule-&gt;n_points();
        
          Point centroid = c.elem-&gt;centroid();
          if( (0.7 &lt;= centroid(0)) && (centroid(0) &lt;= 0.8) &&
              (0.2 &lt;= centroid(1)) && (centroid(1) &lt;= 0.3) )
            for (unsigned int qp=0; qp != n_qpoints; qp++)
              for (unsigned int i=0; i != n_u_dofs; i++)
                c.elem_residual(i) += JxW[qp] * ( 1.*phi[i][qp] ) / 0.01;
        }
        
</pre>
</div>
<div class = "comment">
Build a list (element-by-element) of the Dirichlet dofs. In this case
all boundary dofs are Dirichlet.
</div>

<div class ="fragment">
<pre>
        void initialize_dirichlet_dofs(FEMContext &c, System& system,
                                       std::set&lt;unsigned int&gt;& dirichlet_dofs_set)
        {
          SimpleRB& rb_system = libmesh_cast_ref&lt;SimpleRB&&gt;(system);
          const unsigned int u_var = rb_system.u_var;
        
          std::vector&lt;unsigned int&gt; side_dofs;
          FEInterface::dofs_on_side(c.elem, c.dim, c.element_fe_var[u_var]-&gt;get_fe_type(),
                                    c.side, side_dofs);
        
          for(unsigned int ii=0; ii&lt;side_dofs.size(); ii++)
            dirichlet_dofs_set.insert(c.dof_indices[side_dofs[ii]]);
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
  #include <B><FONT COLOR="#BC8F8F">&quot;gmv_io.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;equation_systems.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;fe.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;quadrature_gauss.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;dof_map.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;sparse_matrix.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;numeric_vector.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;dense_matrix.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;dense_vector.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;fe_interface.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;getpot.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;o_string_stream.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;elem.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;simple_rb.h&quot;</FONT></B>
  
  
  
  
  
  <B><FONT COLOR="#228B22">void</FONT></B> initialize_dirichlet_dofs(FEMContext &amp;c,
                                 System&amp; system,
                                 <B><FONT COLOR="#5F9EA0">std</FONT></B>::set&lt;<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B>&gt;&amp; dirichlet_dofs_set);
  
  Number theta_a_0(RBThetaData&amp; ) { <B><FONT COLOR="#A020F0">return</FONT></B> 0.05; }
  <B><FONT COLOR="#228B22">void</FONT></B> A0(FEMContext&amp;, System&amp;);
  Number theta_a_1(RBThetaData&amp; data) { <B><FONT COLOR="#A020F0">return</FONT></B> data.get_mu()[0]; }
  <B><FONT COLOR="#228B22">void</FONT></B> A1(FEMContext&amp;, System&amp;);
  Number theta_a_2(RBThetaData&amp; data) { <B><FONT COLOR="#A020F0">return</FONT></B> data.get_mu()[1]; }
  <B><FONT COLOR="#228B22">void</FONT></B> A2(FEMContext&amp;, System&amp;);
  
  Number theta_f_0(RBThetaData&amp; ) { <B><FONT COLOR="#A020F0">return</FONT></B> 1.; }
  <B><FONT COLOR="#228B22">void</FONT></B> F0(FEMContext&amp;, System&amp;);
  
  Number theta_l(RBThetaData&amp; ) { <B><FONT COLOR="#A020F0">return</FONT></B> 1.; }
  <B><FONT COLOR="#228B22">void</FONT></B> L0_assembly(FEMContext&amp;, System&amp;);
  <B><FONT COLOR="#228B22">void</FONT></B> L1_assembly(FEMContext&amp;, System&amp;);
  <B><FONT COLOR="#228B22">void</FONT></B> L2_assembly(FEMContext&amp;, System&amp;);
  <B><FONT COLOR="#228B22">void</FONT></B> L3_assembly(FEMContext&amp;, System&amp;);
  
  <B><FONT COLOR="#228B22">int</FONT></B> main (<B><FONT COLOR="#228B22">int</FONT></B> argc, <B><FONT COLOR="#228B22">char</FONT></B>** argv)
  {
    LibMeshInit init (argc, argv);
  
  #<B><FONT COLOR="#A020F0">if</FONT></B> !defined(LIBMESH_HAVE_XDR)
    libmesh_example_assert(false, <B><FONT COLOR="#BC8F8F">&quot;--enable-xdr&quot;</FONT></B>);
  #elif defined(LIBMESH_DEFAULT_SINGLE_PRECISION)
    libmesh_example_assert(false, <B><FONT COLOR="#BC8F8F">&quot;--disable-singleprecision&quot;</FONT></B>);
  #elif !defined(LIBMESH_HAVE_PETSC)
    libmesh_example_assert(false, <B><FONT COLOR="#BC8F8F">&quot;--enable-petsc&quot;</FONT></B>);
  #endif
  
    libmesh_example_assert(2 &lt;= LIBMESH_DIM, <B><FONT COLOR="#BC8F8F">&quot;2D support&quot;</FONT></B>);
    
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::string parameters_filename = <B><FONT COLOR="#BC8F8F">&quot;ex23.in&quot;</FONT></B>;
    GetPot infile(parameters_filename);
  
    <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_elem = infile(<B><FONT COLOR="#BC8F8F">&quot;n_elem&quot;</FONT></B>, 1);       <I><FONT COLOR="#B22222">// Determines the number of elements in the &quot;truth&quot; mesh
</FONT></I>    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> dim = 2;                      <I><FONT COLOR="#B22222">// The number of spatial dimensions
</FONT></I>    
    <B><FONT COLOR="#228B22">bool</FONT></B> online_mode = infile(<B><FONT COLOR="#BC8F8F">&quot;online_mode&quot;</FONT></B>, false); <I><FONT COLOR="#B22222">// Are we in Online mode?
</FONT></I>  
    Mesh mesh (dim);
    <B><FONT COLOR="#5F9EA0">MeshTools</FONT></B>::Generation::build_square (mesh,
                                         n_elem, n_elem,
                                         0., 1.,
                                         0., 1.,
                                         QUAD4);
  
    EquationSystems equation_systems (mesh);
  
    SimpleRB &amp; system =
      equation_systems.add_system&lt;SimpleRB&gt; (<B><FONT COLOR="#BC8F8F">&quot;RBConvectionDiffusion&quot;</FONT></B>);
  
    system.parameters_filename = parameters_filename;
  
    system.attach_dirichlet_dof_initialization(initialize_dirichlet_dofs);
    
    system.attach_A_q(theta_a_0, A0, NULL);
    system.attach_A_q(theta_a_1, A1, NULL);
    system.attach_A_q(theta_a_2, A2, NULL);
    
    system.attach_F_q(theta_f_0, F0, NULL);
    
    system.attach_output(theta_l, L0_assembly, NULL);
    system.attach_output(theta_l, L1_assembly, NULL);
    system.attach_output(theta_l, L2_assembly, NULL);
    system.attach_output(theta_l, L3_assembly, NULL);
    
    system.attach_inner_prod_assembly(A0);
  
    equation_systems.init ();
  
    <B><FONT COLOR="#A020F0">if</FONT></B>(system.initialize_calN_dependent_data)
    {
      mesh.print_info();
      equation_systems.print_info();
    }
    
    system.initialize_RB_system(online_mode);
  
    <B><FONT COLOR="#A020F0">if</FONT></B>(!online_mode)
    {
      system.train_reduced_basis();
      
      system.write_offline_data_to_files();
    }
    <B><FONT COLOR="#A020F0">else</FONT></B>
    {
      system.read_offline_data_from_files();
      
      <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> online_N = infile(<B><FONT COLOR="#BC8F8F">&quot;online_N&quot;</FONT></B>,1);
      <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;Real&gt; online_mu_vector(system.get_n_params());
      <B><FONT COLOR="#A020F0">for</FONT></B>(<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i&lt;system.get_n_params(); i++)
      {
        online_mu_vector[i] = infile(<B><FONT COLOR="#BC8F8F">&quot;online_mu&quot;</FONT></B>, online_mu_vector[i], i);
      }
  
      system.set_current_parameters(online_mu_vector);
      system.print_current_parameters();
  
      system.RB_solve(online_N);
  
      <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;output 1, value = &quot;</FONT></B> &lt;&lt; system.RB_outputs[0]
                &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;, bound = &quot;</FONT></B> &lt;&lt; system.RB_output_error_bounds[0]
                &lt;&lt; std::endl;
      <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;output 2, value = &quot;</FONT></B> &lt;&lt; system.RB_outputs[1]
                &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;, bound = &quot;</FONT></B> &lt;&lt; system.RB_output_error_bounds[1]
                &lt;&lt; std::endl;
      <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;output 3, value = &quot;</FONT></B> &lt;&lt; system.RB_outputs[2]
                &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;, bound = &quot;</FONT></B> &lt;&lt; system.RB_output_error_bounds[2]
                &lt;&lt; std::endl;
      <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;output 4, value = &quot;</FONT></B> &lt;&lt; system.RB_outputs[3]
                &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;, bound = &quot;</FONT></B> &lt;&lt; system.RB_output_error_bounds[3]
                &lt;&lt; std::endl;
  
      <B><FONT COLOR="#A020F0">if</FONT></B>(system.store_basis_functions)
      {
        system.load_RB_solution();
        GMVIO(mesh).write_equation_systems (<B><FONT COLOR="#BC8F8F">&quot;RB_sol.gmv&quot;</FONT></B>,equation_systems);
        
        system.load_basis_function(0);
        GMVIO(mesh).write_equation_systems (<B><FONT COLOR="#BC8F8F">&quot;bf0.gmv&quot;</FONT></B>,equation_systems);
      }
    }
  
    <B><FONT COLOR="#A020F0">return</FONT></B> 0;
  }
  
  <B><FONT COLOR="#228B22">void</FONT></B> A0(FEMContext &amp;c, System&amp; system)
  {
    SimpleRB&amp; rb_system = libmesh_cast_ref&lt;SimpleRB&amp;&gt;(system);
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> u_var = rb_system.u_var;
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Real&gt; &amp;JxW =
      c.element_fe_var[u_var]-&gt;get_JxW();
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;RealGradient&gt; &gt;&amp; dphi =
      c.element_fe_var[u_var]-&gt;get_dphi();
  
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_u_dofs = c.dof_indices_var[u_var].size();
  
    <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_qpoints = c.element_qrule-&gt;n_points();
  
    <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> qp=0; qp != n_qpoints; qp++)
      <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i != n_u_dofs; i++)
        <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j=0; j != n_u_dofs; j++)
          c.elem_jacobian(i,j) += JxW[qp] * dphi[j][qp]*dphi[i][qp];
  }
  
  <B><FONT COLOR="#228B22">void</FONT></B> A1(FEMContext &amp;c, System&amp; system)
  {
    SimpleRB&amp; rb_system = libmesh_cast_ref&lt;SimpleRB&amp;&gt;(system);
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> u_var = rb_system.u_var;
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Real&gt; &amp;JxW =
      c.element_fe_var[u_var]-&gt;get_JxW();
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;Real&gt; &gt;&amp; phi =
      c.element_fe_var[u_var]-&gt;get_phi();
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;RealGradient&gt; &gt;&amp; dphi =
      c.element_fe_var[u_var]-&gt;get_dphi();
  
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_u_dofs = c.dof_indices_var[u_var].size();
  
    <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_qpoints = c.element_qrule-&gt;n_points();
  
    <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> qp=0; qp != n_qpoints; qp++)
      <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i != n_u_dofs; i++)
        <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j=0; j != n_u_dofs; j++)
          c.elem_jacobian(i,j) += JxW[qp] *dphi[j][qp](0)*phi[i][qp];
  }
  
  <B><FONT COLOR="#228B22">void</FONT></B> A2(FEMContext &amp;c, System&amp; system)
  {
    SimpleRB&amp; rb_system = libmesh_cast_ref&lt;SimpleRB&amp;&gt;(system);
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> u_var = rb_system.u_var;
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Real&gt; &amp;JxW =
      c.element_fe_var[u_var]-&gt;get_JxW();
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;Real&gt; &gt;&amp; phi =
      c.element_fe_var[u_var]-&gt;get_phi();
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;RealGradient&gt; &gt;&amp; dphi =
      c.element_fe_var[u_var]-&gt;get_dphi();
  
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_u_dofs = c.dof_indices_var[u_var].size();
  
    <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_qpoints = c.element_qrule-&gt;n_points();
  
    <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> qp=0; qp != n_qpoints; qp++)
      <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i != n_u_dofs; i++)
        <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j=0; j != n_u_dofs; j++)
          c.elem_jacobian(i,j) += JxW[qp] *dphi[j][qp](1)*phi[i][qp];
  }
  
  <B><FONT COLOR="#228B22">void</FONT></B> F0(FEMContext &amp;c, System&amp; system)
  {
    SimpleRB&amp; rb_system = libmesh_cast_ref&lt;SimpleRB&amp;&gt;(system);
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> u_var = rb_system.u_var;
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Real&gt; &amp;JxW =
      c.element_fe_var[u_var]-&gt;get_JxW();
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;Real&gt; &gt;&amp; phi =
      c.element_fe_var[u_var]-&gt;get_phi();
  
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_u_dofs = c.dof_indices_var[u_var].size();
  
    <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_qpoints = c.element_qrule-&gt;n_points();
  
    <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> qp=0; qp != n_qpoints; qp++)
      <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i != n_u_dofs; i++)
        c.elem_residual(i) += JxW[qp] * ( 1.*phi[i][qp] );
  }
  
  <B><FONT COLOR="#228B22">void</FONT></B> L0_assembly(FEMContext &amp;c, System&amp; system)
  {
    SimpleRB&amp; rb_system = libmesh_cast_ref&lt;SimpleRB&amp;&gt;(system);
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> u_var = rb_system.u_var;
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Real&gt; &amp;JxW =
      c.element_fe_var[u_var]-&gt;get_JxW();
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;Real&gt; &gt;&amp; phi =
      c.element_fe_var[u_var]-&gt;get_phi();
  
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_u_dofs = c.dof_indices_var[u_var].size();
  
    <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_qpoints = c.element_qrule-&gt;n_points();
  
    Point centroid = c.elem-&gt;centroid();
    <B><FONT COLOR="#A020F0">if</FONT></B>( (0.7 &lt;= centroid(0)) &amp;&amp; (centroid(0) &lt;= 0.8) &amp;&amp;
        (0.7 &lt;= centroid(1)) &amp;&amp; (centroid(1) &lt;= 0.8) )
      <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> qp=0; qp != n_qpoints; qp++)
        <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i != n_u_dofs; i++)
          c.elem_residual(i) += JxW[qp] * ( 1.*phi[i][qp] ) / 0.01;
  }
  
  <B><FONT COLOR="#228B22">void</FONT></B> L1_assembly(FEMContext &amp;c, System&amp; system)
  {
    SimpleRB&amp; rb_system = libmesh_cast_ref&lt;SimpleRB&amp;&gt;(system);
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> u_var = rb_system.u_var;
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Real&gt; &amp;JxW =
      c.element_fe_var[u_var]-&gt;get_JxW();
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;Real&gt; &gt;&amp; phi =
      c.element_fe_var[u_var]-&gt;get_phi();
  
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_u_dofs = c.dof_indices_var[u_var].size();
  
    <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_qpoints = c.element_qrule-&gt;n_points();
  
    Point centroid = c.elem-&gt;centroid();
    <B><FONT COLOR="#A020F0">if</FONT></B>( (0.2 &lt;= centroid(0)) &amp;&amp; (centroid(0) &lt;= 0.3) &amp;&amp;
        (0.7 &lt;= centroid(1)) &amp;&amp; (centroid(1) &lt;= 0.8) )
      <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> qp=0; qp != n_qpoints; qp++)
        <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i != n_u_dofs; i++)
          c.elem_residual(i) += JxW[qp] * ( 1.*phi[i][qp] ) / 0.01;
  }
  
  <B><FONT COLOR="#228B22">void</FONT></B> L2_assembly(FEMContext &amp;c, System&amp; system)
  {
    SimpleRB&amp; rb_system = libmesh_cast_ref&lt;SimpleRB&amp;&gt;(system);
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> u_var = rb_system.u_var;
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Real&gt; &amp;JxW =
      c.element_fe_var[u_var]-&gt;get_JxW();
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;Real&gt; &gt;&amp; phi =
      c.element_fe_var[u_var]-&gt;get_phi();
  
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_u_dofs = c.dof_indices_var[u_var].size();
  
    <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_qpoints = c.element_qrule-&gt;n_points();
  
    Point centroid = c.elem-&gt;centroid();
    <B><FONT COLOR="#A020F0">if</FONT></B>( (0.2 &lt;= centroid(0)) &amp;&amp; (centroid(0) &lt;= 0.3) &amp;&amp;
        (0.2 &lt;= centroid(1)) &amp;&amp; (centroid(1) &lt;= 0.3) )
      <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> qp=0; qp != n_qpoints; qp++)
        <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i != n_u_dofs; i++)
          c.elem_residual(i) += JxW[qp] * ( 1.*phi[i][qp] ) / 0.01;
  }
  
  <B><FONT COLOR="#228B22">void</FONT></B> L3_assembly(FEMContext &amp;c, System&amp; system)
  {
    SimpleRB&amp; rb_system = libmesh_cast_ref&lt;SimpleRB&amp;&gt;(system);
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> u_var = rb_system.u_var;
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Real&gt; &amp;JxW =
      c.element_fe_var[u_var]-&gt;get_JxW();
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;Real&gt; &gt;&amp; phi =
      c.element_fe_var[u_var]-&gt;get_phi();
  
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_u_dofs = c.dof_indices_var[u_var].size();
  
    <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_qpoints = c.element_qrule-&gt;n_points();
  
    Point centroid = c.elem-&gt;centroid();
    <B><FONT COLOR="#A020F0">if</FONT></B>( (0.7 &lt;= centroid(0)) &amp;&amp; (centroid(0) &lt;= 0.8) &amp;&amp;
        (0.2 &lt;= centroid(1)) &amp;&amp; (centroid(1) &lt;= 0.3) )
      <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> qp=0; qp != n_qpoints; qp++)
        <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i != n_u_dofs; i++)
          c.elem_residual(i) += JxW[qp] * ( 1.*phi[i][qp] ) / 0.01;
  }
  
  <B><FONT COLOR="#228B22">void</FONT></B> initialize_dirichlet_dofs(FEMContext &amp;c, System&amp; system,
                                 <B><FONT COLOR="#5F9EA0">std</FONT></B>::set&lt;<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B>&gt;&amp; dirichlet_dofs_set)
  {
    SimpleRB&amp; rb_system = libmesh_cast_ref&lt;SimpleRB&amp;&gt;(system);
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> u_var = rb_system.u_var;
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B>&gt; side_dofs;
    <B><FONT COLOR="#5F9EA0">FEInterface</FONT></B>::dofs_on_side(c.elem, c.dim, c.element_fe_var[u_var]-&gt;get_fe_type(),
                              c.side, side_dofs);
  
    <B><FONT COLOR="#A020F0">for</FONT></B>(<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> ii=0; ii&lt;side_dofs.size(); ii++)
      dirichlet_dofs_set.insert(c.dof_indices[side_dofs[ii]]);
  }
</pre> 
<a name="output"></a> 
<br><br><br> <h1> The console output of the program: </h1> 
<pre>
Updated .depend
***************************************************************
* Running  ./ex23-opt
***************************************************************
 
*** Warning, This code is untested, experimental, or likely to see future API changes: src/systems/rb_base.C, line 68, compiled Feb  1 2011 at 13:24:01 ***

RBSystem parameters:
system name: RBConvectionDiffusion
constrained_problem: 0
Nmax: 20
A_q operators attached: 3
F_q functions attached: 1
Number of A EIM systems: 0
Number of F EIM systems: 0
n_outputs: 4
output 0, Q_l = 1
output 1, Q_l = 1
output 2, Q_l = 1
output 3, Q_l = 1
Parameter 0: Min = -2, Max = 2, log scaling = 0
Parameter 1: Min = -2, Max = 2, log scaling = 0
n_training_samples: 100
using deterministic training samples? 1
store/load basis functions? 1
  write out basis functions in binary format? 1
  read in basis functions in binary format? 1
store/load residual representors? 0
low-memory mode? 0
reuse preconditioner? 1
return a relative error bound from RB_solve? 1
write out data during basis training? 0
initializing calN-dependent data structures? 1
impose internal Dirichlet BCs? 0
impose internal fluxes? 0
quiet mode? 1
initial parameter: mu[0] = 2, mu[1] = 2

 Mesh Information:
  mesh_dimension()=2
  spatial_dimension()=3
  n_nodes()=2601
    n_local_nodes()=2601
  n_elem()=2500
    n_local_elem()=2500
    n_active_elem()=2500
  n_subdomains()=1
  n_processors()=1
  processor_id()=0

 EquationSystems
  n_systems()=1
   System "RBConvectionDiffusion"
    Type "RBSystem"
    Variables="u" 
    Finite Element Types="LAGRANGE", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" 
    Approximation Orders="FIRST", "THIRD" 
    n_dofs()=2601
    n_local_dofs()=2601
    n_constrained_dofs()=0
    n_vectors()=1

Compute output dual norms
output_dual_norms[0][0] = 0.323675
output_dual_norms[1][0] = 0.323675
output_dual_norms[2][0] = 0.323675
output_dual_norms[3][0] = 0.323675

---- Training solve 1 ----
mu[0] = 2
mu[1] = 2


Enriching the RB space
Reduced basis dimension = 1
Updating RB matrices
Updating RB residual terms
Performing RB solves on training set
Maximum a posteriori error is 0.875694


---- Training solve 2 ----
mu[0] = -2
mu[1] = -2


Enriching the RB space
Reduced basis dimension = 2
Updating RB matrices
Updating RB residual terms
Performing RB solves on training set
Maximum a posteriori error is 0.510877


---- Training solve 3 ----
mu[0] = -2
mu[1] = 2


Enriching the RB space
Reduced basis dimension = 3
Updating RB matrices
Updating RB residual terms
Performing RB solves on training set
Maximum a posteriori error is 0.47092


---- Training solve 4 ----
mu[0] = 2
mu[1] = -2


Enriching the RB space
Reduced basis dimension = 4
Updating RB matrices
Updating RB residual terms
Performing RB solves on training set
Maximum a posteriori error is 0.245027


---- Training solve 5 ----
mu[0] = 0.222222
mu[1] = 0.222222


Enriching the RB space
Reduced basis dimension = 5
Updating RB matrices
Updating RB residual terms
Performing RB solves on training set
Maximum a posteriori error is 0.193299


---- Training solve 6 ----
mu[0] = -0.666667
mu[1] = -0.666667


Enriching the RB space
Reduced basis dimension = 6
Updating RB matrices
Updating RB residual terms
Performing RB solves on training set
Maximum a posteriori error is 0.176492


---- Training solve 7 ----
mu[0] = -0.222222
mu[1] = 2


Enriching the RB space
Reduced basis dimension = 7
Updating RB matrices
Updating RB residual terms
Performing RB solves on training set
Maximum a posteriori error is 0.173431


---- Training solve 8 ----
mu[0] = 2
mu[1] = -0.222222


Enriching the RB space
Reduced basis dimension = 8
Updating RB matrices
Updating RB residual terms
Performing RB solves on training set
Maximum a posteriori error is 0.146806


---- Training solve 9 ----
mu[0] = -2
mu[1] = 0.222222


Enriching the RB space
Reduced basis dimension = 9
Updating RB matrices
Updating RB residual terms
Performing RB solves on training set
Maximum a posteriori error is 0.13704


---- Training solve 10 ----
mu[0] = 0.666667
mu[1] = -1.11111


Enriching the RB space
Reduced basis dimension = 10
Updating RB matrices
Updating RB residual terms
Performing RB solves on training set
Maximum a posteriori error is 0.10405


---- Training solve 11 ----
mu[0] = -0.222222
mu[1] = -2


Enriching the RB space
Reduced basis dimension = 11
Updating RB matrices
Updating RB residual terms
Performing RB solves on training set
Maximum a posteriori error is 0.0896627


---- Training solve 12 ----
mu[0] = -0.666667
mu[1] = 0.666667


Enriching the RB space
Reduced basis dimension = 12
Updating RB matrices
Updating RB residual terms
Performing RB solves on training set
Maximum a posteriori error is 0.0612686


---- Training solve 13 ----
mu[0] = 0.666667
mu[1] = 1.11111


Enriching the RB space
Reduced basis dimension = 13
Updating RB matrices
Updating RB residual terms
Performing RB solves on training set
Maximum a posteriori error is 0.0488624


---- Training solve 14 ----
mu[0] = 0.666667
mu[1] = -0.222222


Enriching the RB space
Reduced basis dimension = 14
Updating RB matrices
Updating RB residual terms
Performing RB solves on training set
Maximum a posteriori error is 0.0459427


---- Training solve 15 ----
mu[0] = -0.222222
mu[1] = -0.222222


Enriching the RB space
Reduced basis dimension = 15
Updating RB matrices
Updating RB residual terms
Performing RB solves on training set
Maximum a posteriori error is 0.0325616


---- Training solve 16 ----
mu[0] = -1.11111
mu[1] = 1.55556


Enriching the RB space
Reduced basis dimension = 16
Updating RB matrices
Updating RB residual terms
Performing RB solves on training set
Maximum a posteriori error is 0.0315341


---- Training solve 17 ----
mu[0] = 1.11111
mu[1] = 0.666667


Enriching the RB space
Reduced basis dimension = 17
Updating RB matrices
Updating RB residual terms
Performing RB solves on training set
Maximum a posteriori error is 0.0291748


---- Training solve 18 ----
mu[0] = 0.222222
mu[1] = -0.666667


Enriching the RB space
Reduced basis dimension = 18
Updating RB matrices
Updating RB residual terms
Performing RB solves on training set
Maximum a posteriori error is 0.0273367


---- Training solve 19 ----
mu[0] = -2
mu[1] = -0.666667


Enriching the RB space
Reduced basis dimension = 19
Updating RB matrices
Updating RB residual terms
Performing RB solves on training set
Maximum a posteriori error is 0.0266323


---- Training solve 20 ----
mu[0] = -1.11111
mu[1] = 0.222222


Enriching the RB space
Reduced basis dimension = 20
Updating RB matrices
Updating RB residual terms
Performing RB solves on training set
Maximum a posteriori error is 0.0222554

Maximum number of basis functions reached: Nmax = 20
In RBSystem::write_offline_data_to_files, directory offline_data already exists, overwriting contents.
Writing out the basis functions...

-------------------------------------------------------------------
| Time:           Thu Feb  3 12:13:01 2011                         |
| OS:             Linux                                            |
| HostName:       daedalus                                         |
| OS Release:     2.6.32-26-generic                                |
| OS Version:     #46-Ubuntu SMP Tue Oct 26 16:47:18 UTC 2010      |
| Machine:        x86_64                                           |
| Username:       roystgnr                                         |
| Configuration:  ./configure run on Tue Feb  1 12:58:27 CST 2011  |
-------------------------------------------------------------------
 ---------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=1.08326, Active time=1.06499                                                  |
 ---------------------------------------------------------------------------------------------------------------
| Event                             nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                             w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|---------------------------------------------------------------------------------------------------------------|
|                                                                                                               |
|                                                                                                               |
| DofMap                                                                                                        |
|   add_neighbors_to_send_list()    1         0.0005      0.000452    0.0005      0.000452    0.04     0.04     |
|   compute_sparsity()              1         0.0030      0.002950    0.0038      0.003843    0.28     0.36     |
|   create_dof_constraints()        1         0.0003      0.000272    0.0003      0.000272    0.03     0.03     |
|   distribute_dofs()               1         0.0009      0.000902    0.0028      0.002773    0.08     0.26     |
|   dof_indices()                   52500     0.0150      0.000000    0.0150      0.000000    1.41     1.41     |
|   prepare_send_list()             1         0.0000      0.000000    0.0000      0.000000    0.00     0.00     |
|   reinit()                        1         0.0019      0.001869    0.0019      0.001869    0.18     0.18     |
|                                                                                                               |
| FE                                                                                                            |
|   compute_affine_map()            25200     0.0123      0.000000    0.0123      0.000000    1.16     1.16     |
|   compute_face_map()              200       0.0004      0.000002    0.0008      0.000004    0.04     0.08     |
|   compute_shape_functions()       25200     0.0060      0.000000    0.0060      0.000000    0.56     0.56     |
|   init_face_shape_functions()     105       0.0001      0.000000    0.0001      0.000000    0.00     0.00     |
|   init_shape_functions()          210       0.0007      0.000003    0.0007      0.000003    0.07     0.07     |
|   inverse_map()                   800       0.0007      0.000001    0.0007      0.000001    0.07     0.07     |
|                                                                                                               |
| Mesh                                                                                                          |
|   find_neighbors()                1         0.0030      0.002984    0.0030      0.002984    0.28     0.28     |
|   renumber_nodes_and_elem()       2         0.0001      0.000056    0.0001      0.000056    0.01     0.01     |
|                                                                                                               |
| MeshTools::Generation                                                                                         |
|   build_cube()                    1         0.0007      0.000746    0.0007      0.000746    0.07     0.07     |
|                                                                                                               |
| Parallel                                                                                                      |
|   allgather()                     1         0.0000      0.000001    0.0000      0.000001    0.00     0.00     |
|   receive()                       80        0.0001      0.000001    0.0001      0.000001    0.01     0.01     |
|   send()                          80        0.0002      0.000002    0.0002      0.000002    0.01     0.01     |
|                                                                                                               |
| Partitioner                                                                                                   |
|   single_partition()              1         0.0003      0.000253    0.0003      0.000253    0.02     0.02     |
|                                                                                                               |
| PetscLinearSolver                                                                                             |
|   solve()                         85        0.5664      0.006663    0.5664      0.006663    53.18    53.18    |
|                                                                                                               |
| RBSystem                                                                                                      |
|   RB_solve()                      2000      0.0174      0.000009    0.0752      0.000038    1.64     7.06     |
|   add_scaled_matrix_and_vector()  9         0.0625      0.006940    0.0978      0.010863    5.87     9.18     |
|   compute_max_error_bound()       20        0.0019      0.000093    0.0772      0.003858    0.17     7.24     |
|   compute_output_dual_norms()     1         0.0047      0.004748    0.0357      0.035689    0.45     3.35     |
|   compute_residual_dual_norm()    2000      0.0576      0.000029    0.0576      0.000029    5.41     5.41     |
|   enrich_RB_space()               20        0.0127      0.000635    0.0127      0.000635    1.19     1.19     |
|   initialize_dirichlet_dofs()     1         0.0072      0.007170    0.0131      0.013078    0.67     1.23     |
|   train_reduced_basis()           1         0.0013      0.001349    0.8949      0.894937    0.13     84.03    |
|   truth_assembly()                20        0.0218      0.001090    0.0218      0.001090    2.05     2.05     |
|   truth_solve()                   20        0.0023      0.000114    0.1177      0.005884    0.21     11.05    |
|   update_RB_system_matrices()     20        0.0651      0.003256    0.0651      0.003256    6.11     6.11     |
|   update_residual_terms()         20        0.1424      0.007118    0.5852      0.029261    13.37    54.95    |
|   write_offline_data_to_files()   1         0.0547      0.054738    0.0550      0.054978    5.14     5.16     |
|   zero_dirichlet_dofs_on_rhs()    65        0.0001      0.000001    0.0010      0.000015    0.01     0.09     |
|   zero_dirichlet_dofs_on_vector() 65        0.0009      0.000014    0.0009      0.000014    0.09     0.09     |
 ---------------------------------------------------------------------------------------------------------------
| Totals:                           108735    1.0650                                          100.00            |
 ---------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running  ./ex23-opt
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
