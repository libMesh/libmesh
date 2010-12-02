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
Compiling C++ (in optimized mode) ex23.C...
/org/centers/pecos/LIBRARIES/GCC/gcc-4.5.1-lucid/libexec/gcc/x86_64-unknown-linux-gnu/4.5.1/cc1plus: error while loading shared libraries: libmpc.so.2: cannot open shared object file: No such file or directory
make[1]: *** [ex23.x86_64-unknown-linux-gnu.opt.o] Error 1
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
