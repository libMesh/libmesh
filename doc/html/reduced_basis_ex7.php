<?php $root=""; ?>
<?php require($root."navigation.php"); ?>
<html>
<head>
  <?php load_style($root); ?>
</head>
 
<body>
 
<?php make_navigation("reduced_basis_ex7",$root)?>
 
<div class="content">
<a name="comments"></a> 
<br><br><br> <h1> The source file assembly.h with comments: </h1> 
<div class = "comment">
</div>

<div class ="fragment">
<pre>
        #ifndef __assembly_h__
        #define __assembly_h__
        
        #include "libmesh/sparse_matrix.h"
        #include "libmesh/numeric_vector.h"
        #include "libmesh/dense_matrix.h"
        #include "libmesh/dense_vector.h"
        #include "libmesh/fe.h"
        #include "libmesh/fe_interface.h"
        #include "libmesh/fe_base.h"
        #include "libmesh/elem_assembly.h"
        #include "libmesh/quadrature_gauss.h"
        
</pre>
</div>
<div class = "comment">
rbOOmit includes
</div>

<div class ="fragment">
<pre>
        #include "libmesh/rb_theta.h"
        #include "libmesh/rb_assembly_expansion.h"
        
        #define damping_epsilon 0.001
        #define R_rad 12.0
        
        #ifdef LIBMESH_USE_COMPLEX_NUMBERS
        
</pre>
</div>
<div class = "comment">
Bring in bits from the libMesh namespace.
Just the bits we're using, since this is a header.
</div>

<div class ="fragment">
<pre>
        using libMesh::ElemAssembly;
        using libMesh::FEInterface;
        using libMesh::FEMContext;
        using libMesh::Number;
        using libMesh::Point;
        using libMesh::RBAssemblyExpansion;
        using libMesh::RBParameters;
        using libMesh::RBTheta;
        using libMesh::RBThetaExpansion;
        using libMesh::Real;
        using libMesh::RealGradient;
        using libMesh::MeshBase;
        using libMesh::libmesh_conj;
        
</pre>
</div>
<div class = "comment">
Functors for the parameter-dependent part of the affine decomposition of the PDE
</div>

<div class ="fragment">
<pre>
        struct ThetaA0 : RBTheta { virtual Number evaluate(const RBParameters& mu) { Number val(1., mu.get_value("frequency")*damping_epsilon); return val; } };
        struct ThetaA1 : RBTheta { virtual Number evaluate(const RBParameters& mu) { Number val(-mu.get_value("frequency")*mu.get_value("frequency"), 0.); return val; } };
        struct ThetaA2 : RBTheta { virtual Number evaluate(const RBParameters& mu) { Number val(0., mu.get_value("frequency")); return val; } };
        struct ThetaA3 : RBTheta { virtual Number evaluate(const RBParameters& mu) { Number val(0.5/R_rad, mu.get_value("frequency")); return val; } };
        
        struct ThetaF0 : RBTheta { virtual Number evaluate(const RBParameters& mu) { Number val(0., 2.*mu.get_value("frequency")); return val; } };
        
        struct ThetaOutput0 : RBTheta { virtual Number evaluate(const RBParameters& ) { Number val(1., 0.); return val; } };
        
        struct AcousticsInnerProduct : ElemAssembly
        {
          virtual void interior_assembly(FEMContext &c)
          {
            const unsigned int p_var = 0;
        
            const std::vector&lt;Real&gt; &JxW =
              c.element_fe_var[p_var]-&gt;get_JxW();
        
            const std::vector&lt;std::vector&lt;Real&gt; &gt;& phi =
              c.element_fe_var[p_var]-&gt;get_phi();
              
            const std::vector&lt;std::vector&lt;RealGradient&gt; &gt;& dphi =
              c.element_fe_var[p_var]-&gt;get_dphi();
        
</pre>
</div>
<div class = "comment">
The number of local degrees of freedom in each variable
</div>

<div class ="fragment">
<pre>
            const unsigned int n_p_dofs = c.dof_indices_var[p_var].size();
        
</pre>
</div>
<div class = "comment">
Now we will build the affine operator
</div>

<div class ="fragment">
<pre>
            unsigned int n_qpoints = (c.get_element_qrule())-&gt;n_points();
        
            for (unsigned int qp=0; qp != n_qpoints; qp++)
              for (unsigned int i=0; i != n_p_dofs; i++)
                for (unsigned int j=0; j != n_p_dofs; j++)
                  c.elem_jacobian(i,j) += JxW[qp] * (dphi[j][qp](0)*libmesh_conj(dphi[i][qp](0)) +
                                                     dphi[j][qp](1)*libmesh_conj(dphi[i][qp](1)) + 
                                                     (phi[j][qp]*libmesh_conj(phi[i][qp])) );
          }
        };
        
        struct A0 : ElemAssembly
        {
          virtual void interior_assembly(FEMContext &c)
          {
            const unsigned int p_var = 0;
        
            const std::vector&lt;Real&gt; &JxW =
              c.element_fe_var[p_var]-&gt;get_JxW();
        
</pre>
</div>
<div class = "comment">
The velocity shape function gradients at interior
quadrature points.
</div>

<div class ="fragment">
<pre>
            const std::vector&lt;std::vector&lt;RealGradient&gt; &gt;& dphi =
              c.element_fe_var[p_var]-&gt;get_dphi();
        
</pre>
</div>
<div class = "comment">
The number of local degrees of freedom in each variable
</div>

<div class ="fragment">
<pre>
            const unsigned int n_p_dofs = c.dof_indices_var[p_var].size();
        
</pre>
</div>
<div class = "comment">
Now we will build the affine operator
</div>

<div class ="fragment">
<pre>
            unsigned int n_qpoints = (c.get_element_qrule())-&gt;n_points();
        
            for (unsigned int qp=0; qp != n_qpoints; qp++)
              for (unsigned int i=0; i != n_p_dofs; i++)
                for (unsigned int j=0; j != n_p_dofs; j++)
                  c.elem_jacobian(i,j) += JxW[qp] * (dphi[j][qp](0)*libmesh_conj(dphi[i][qp](0)) +
                                                     dphi[j][qp](1)*libmesh_conj(dphi[i][qp](1)));
          }
        };
        
        
        struct A1 : ElemAssembly
        {
          virtual void interior_assembly(FEMContext &c)
          {
            const unsigned int p_var = 0;
        
            const std::vector&lt;Real&gt; &JxW =
              c.element_fe_var[p_var]-&gt;get_JxW();
        
            const std::vector&lt;std::vector&lt;Real&gt; &gt;& phi =
              c.element_fe_var[p_var]-&gt;get_phi();
        
</pre>
</div>
<div class = "comment">
The number of local degrees of freedom in each variable
</div>

<div class ="fragment">
<pre>
            const unsigned int n_p_dofs = c.dof_indices_var[p_var].size();
        
</pre>
</div>
<div class = "comment">
Now we will build the affine operator
</div>

<div class ="fragment">
<pre>
            unsigned int n_qpoints = (c.get_element_qrule())-&gt;n_points();
        
            for (unsigned int qp=0; qp != n_qpoints; qp++)
              for (unsigned int i=0; i != n_p_dofs; i++)
                for (unsigned int j=0; j != n_p_dofs; j++)
                  c.elem_jacobian(i,j) += JxW[qp] * (phi[j][qp]*libmesh_conj(phi[i][qp]));
          }
        };
        
        struct A2 : ElemAssembly
        {
          virtual void boundary_assembly(FEMContext &c)
          {
            if( mesh-&gt;boundary_info-&gt;has_boundary_id (c.elem, c.side, 1) ) // Forcing on the horn "inlet"
            {
              const unsigned int p_var = 0;
        
              const std::vector&lt;Real&gt; &JxW_face =
                c.side_fe_var[p_var]-&gt;get_JxW();
        
              const std::vector&lt;std::vector&lt;Real&gt; &gt;& phi_face =
                c.side_fe_var[p_var]-&gt;get_phi();
        
</pre>
</div>
<div class = "comment">
The number of local degrees of freedom in each variable
</div>

<div class ="fragment">
<pre>
              const unsigned int n_p_dofs = c.dof_indices_var[p_var].size();
        
</pre>
</div>
<div class = "comment">
Now we will build the affine operator
</div>

<div class ="fragment">
<pre>
              unsigned int n_sidepoints = (c.side_qrule)-&gt;n_points();
        
              for (unsigned int qp=0; qp != n_sidepoints; qp++)
                for (unsigned int i=0; i != n_p_dofs; i++)
                  for (unsigned int j=0; j != n_p_dofs; j++)
                    c.elem_jacobian(i,j) += JxW_face[qp] * phi_face[j][qp] * libmesh_conj(phi_face[i][qp]);
            }
          }
          
          MeshBase* mesh;
        };
        
        struct A3 : ElemAssembly
        {
          virtual void boundary_assembly(FEMContext &c)
          {
            if( mesh-&gt;boundary_info-&gt;has_boundary_id (c.elem, c.side, 2) ) // Radiation condition on the "bubble"
            {
              const unsigned int p_var = 0;
        
              const std::vector&lt;Real&gt; &JxW_face =
                c.side_fe_var[p_var]-&gt;get_JxW();
        
              const std::vector&lt;std::vector&lt;Real&gt; &gt;& phi_face =
                c.side_fe_var[p_var]-&gt;get_phi();
        
</pre>
</div>
<div class = "comment">
The number of local degrees of freedom in each variable
</div>

<div class ="fragment">
<pre>
              const unsigned int n_p_dofs = c.dof_indices_var[p_var].size();
        
</pre>
</div>
<div class = "comment">
Now we will build the affine operator
</div>

<div class ="fragment">
<pre>
              unsigned int n_sidepoints = (c.side_qrule)-&gt;n_points();
        
              for (unsigned int qp=0; qp != n_sidepoints; qp++)
                for (unsigned int i=0; i != n_p_dofs; i++)
                  for (unsigned int j=0; j != n_p_dofs; j++)
                    c.elem_jacobian(i,j) += JxW_face[qp] * phi_face[j][qp] * libmesh_conj(phi_face[i][qp]);
            }
          }
          
          MeshBase* mesh;
        };
        
        struct F0 : ElemAssembly
        {
          virtual void boundary_assembly(FEMContext &c)
          {
            if( mesh-&gt;boundary_info-&gt;has_boundary_id (c.elem, c.side, 1) ) // Output is calculated on the horn "inlet"
            {
              const unsigned int p_var = 0;
        
              const std::vector&lt;Real&gt; &JxW_face =
                c.side_fe_var[p_var]-&gt;get_JxW();
        
              const std::vector&lt;std::vector&lt;Real&gt; &gt;& phi_face =
                c.side_fe_var[p_var]-&gt;get_phi();
        
</pre>
</div>
<div class = "comment">
The number of local degrees of freedom in each variable
</div>

<div class ="fragment">
<pre>
              const unsigned int n_p_dofs = c.dof_indices_var[p_var].size();
        
</pre>
</div>
<div class = "comment">
Now we will build the affine operator
</div>

<div class ="fragment">
<pre>
              unsigned int n_sidepoints = (c.side_qrule)-&gt;n_points();
        
              for (unsigned int qp=0; qp != n_sidepoints; qp++)
                for (unsigned int i=0; i != n_p_dofs; i++)
                    c.elem_residual(i) += JxW_face[qp] * libmesh_conj(phi_face[i][qp]);
            }
          }
          
          MeshBase* mesh;
        };
        
        struct Output0 : ElemAssembly
        {
          virtual void boundary_assembly(FEMContext &c)
          {
            if( mesh-&gt;boundary_info-&gt;has_boundary_id (c.elem, c.side, 1) ) // Forcing on the horn "inlet"
            {
              const unsigned int p_var = 0;
        
              const std::vector&lt;Real&gt; &JxW_face =
                c.side_fe_var[p_var]-&gt;get_JxW();
        
              const std::vector&lt;std::vector&lt;Real&gt; &gt;& phi_face =
                c.side_fe_var[p_var]-&gt;get_phi();
        
</pre>
</div>
<div class = "comment">
The number of local degrees of freedom in each variable
</div>

<div class ="fragment">
<pre>
              const unsigned int n_p_dofs = c.dof_indices_var[p_var].size();
        
</pre>
</div>
<div class = "comment">
Now we will build the affine operator
</div>

<div class ="fragment">
<pre>
              unsigned int n_sidepoints = (c.side_qrule)-&gt;n_points();
        
              for (unsigned int qp=0; qp != n_sidepoints; qp++)
                for (unsigned int i=0; i != n_p_dofs; i++)
                    c.elem_residual(i) += JxW_face[qp] * libmesh_conj(phi_face[i][qp]);
            }
          }
          
          MeshBase* mesh;
        };
        
        struct AcousticsRBThetaExpansion : RBThetaExpansion
        {
        
          /**
           * Constructor.
           */
          AcousticsRBThetaExpansion()
          {
</pre>
</div>
<div class = "comment">
set up the RBThetaExpansion object
</div>

<div class ="fragment">
<pre>
            attach_A_theta(&theta_a_0);   // Attach the lhs theta
            attach_A_theta(&theta_a_1);
            attach_A_theta(&theta_a_2);
            attach_A_theta(&theta_a_3);
        
            attach_F_theta(&theta_f_0);    // Attach the rhs theta
            
            attach_output_theta(&theta_output_0);
          }
        
</pre>
</div>
<div class = "comment">
The RBTheta member variables
</div>

<div class ="fragment">
<pre>
          ThetaA0 theta_a_0;
          ThetaA1 theta_a_1;
          ThetaA2 theta_a_2;
          ThetaA3 theta_a_3;
          ThetaF0 theta_f_0;
          ThetaOutput0 theta_output_0;
        };
        
</pre>
</div>
<div class = "comment">
Define an RBAssemblyExpansion class for this PDE
</div>

<div class ="fragment">
<pre>
        struct AcousticsRBAssemblyExpansion : RBAssemblyExpansion
        {
        
          /**
           * Constructor.
           */
          AcousticsRBAssemblyExpansion(MeshBase& mesh_in)
          {
            A2_assembly.mesh = &mesh_in;
            A3_assembly.mesh = &mesh_in;
            F0_assembly.mesh = &mesh_in;
            Output0_assembly.mesh = &mesh_in;
            
</pre>
</div>
<div class = "comment">
And set up the RBAssemblyExpansion object
</div>

<div class ="fragment">
<pre>
            attach_A_assembly(&A0_assembly); // Attach the lhs assembly
            attach_A_assembly(&A1_assembly);
            attach_A_assembly(&A2_assembly);
            attach_A_assembly(&A3_assembly);
            
            attach_F_assembly(&F0_assembly); // Attach the rhs assembly
            
            attach_output_assembly(&Output0_assembly);
          }
        
</pre>
</div>
<div class = "comment">
The ElemAssembly objects
</div>

<div class ="fragment">
<pre>
          A0 A0_assembly;
          A1 A1_assembly;
          A2 A2_assembly;
          A3 A3_assembly;
          F0 F0_assembly;
          Output0 Output0_assembly;
          AcousticsInnerProduct acoustics_inner_product;
        };
        
        #endif
        
        #endif
</pre>
</div>

<a name="comments"></a> 
<br><br><br> <h1> The source file rb_classes.h with comments: </h1> 
<div class = "comment">
  

<br><br>rbOOmit is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.
  

<br><br>You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


<br><br></div>

<div class ="fragment">
<pre>
        #ifndef __rb_classes_h__
        #define __rb_classes_h__
        
        #include "libmesh/rb_construction.h"
        #include "libmesh/fe_base.h"
        
</pre>
</div>
<div class = "comment">
local include
</div>

<div class ="fragment">
<pre>
        #include "assembly.h"
        
        #ifdef LIBMESH_USE_COMPLEX_NUMBERS
        
</pre>
</div>
<div class = "comment">
Bring in bits from the libMesh namespace.
Just the bits we're using, since this is a header.
</div>

<div class ="fragment">
<pre>
        using libMesh::AutoPtr;
        using libMesh::DirichletBoundary;
        using libMesh::EquationSystems;
        using libMesh::FEMContext;
        using libMesh::RBConstruction;
        using libMesh::RBEvaluation;
        using libMesh::Real;
        
        
</pre>
</div>
<div class = "comment">
A simple subclass of RBEvaluation, which just needs to specify
(a lower bound for) the coercivity constant for this problem.
For this simple convection-diffusion problem, we can set the
coercivity constant lower bound to 0.05.
</div>

<div class ="fragment">
<pre>
        class SimpleRBEvaluation : public RBEvaluation
        {
        public:
        
          /**
           * Constructor. Just set the theta expansion.
           */
          SimpleRBEvaluation()
          {
            set_rb_theta_expansion(acoustics_rb_theta_expansion);
          }
        
          /**
           * Return a "dummy" stability lower bound factor (for a rigorous error bound
           * this should be a lower bound for the frequency dependent inf-sup constant)
           */
          virtual Real get_stability_lower_bound() { return 1.; }
        
          /**
           * The object that stores the "theta" expansion of the parameter dependent PDE,
           * i.e. the set of parameter-dependent functions in the affine expansion of the PDE.
           */
          AcousticsRBThetaExpansion acoustics_rb_theta_expansion;
        
        };
        
</pre>
</div>
<div class = "comment">
A simple subclass of Construction, which just needs to override build_rb_evaluation
in order to build a SimpleRBEvaluation object, rather than an RBEvaluation object.
</div>

<div class ="fragment">
<pre>
        class SimpleRBConstruction : public RBConstruction
        {
        public:
        
          SimpleRBConstruction (EquationSystems& es,
                                const std::string& name,
                                const unsigned int number)
          : Parent(es, name, number)
          {}
        
          /**
           * Destructor.
           */
          virtual ~SimpleRBConstruction () { delete acoustics_rb_assembly_expansion; }
        
          /**
           * The type of system.
           */
          typedef SimpleRBConstruction sys_type;
        
          /**
           * The type of the parent.
           */
          typedef RBConstruction Parent;
        
          /**
           * Initialize data structures.
           */
          virtual void init_data()
          {
            p_var = this-&gt;add_variable ("p", SECOND);
        
            Parent::init_data();
            
            acoustics_rb_assembly_expansion = new AcousticsRBAssemblyExpansion(get_mesh());
        
</pre>
</div>
<div class = "comment">
Set the rb_assembly_expansion for this Construction object.
The theta expansion comes from the RBEvaluation object.
</div>

<div class ="fragment">
<pre>
            set_rb_assembly_expansion(*acoustics_rb_assembly_expansion);
        
</pre>
</div>
<div class = "comment">
We need to define an inner product matrix for this problem
</div>

<div class ="fragment">
<pre>
            set_inner_product_assembly(acoustics_rb_assembly_expansion-&gt;acoustics_inner_product);
          }
        
          /**
           * Pre-request all relevant element data.
           */
          virtual void init_context(FEMContext &c)
          {
</pre>
</div>
<div class = "comment">
For efficiency, we should prerequest all
the data we will need to build the
linear system before doing an element loop.
</div>

<div class ="fragment">
<pre>
            c.element_fe_var[p_var]-&gt;get_JxW();
            c.element_fe_var[p_var]-&gt;get_phi();
            c.element_fe_var[p_var]-&gt;get_dphi();
          }
        
          /**
           * Variable number for pd.
           */
          unsigned int p_var;
          
          /**
           * The object that stores the "assembly" expansion of the parameter dependent PDE,
           * i.e. the objects that define how to assemble the set of parameter-independent
           * operators in the affine expansion of the PDE.
           */
          AcousticsRBAssemblyExpansion* acoustics_rb_assembly_expansion;
        
        };
        
        #endif
        
        #endif
</pre>
</div>

<a name="comments"></a> 
<br><br><br> <h1> The source file reduced_basis_ex7.C with comments: </h1> 
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
<h1>Reduced Basis Example 7 - Acoustic Horn</h1>
 

<br><br>In this example problem we use the Certified Reduced Basis method
to solve a complex-valued Helmholtz problem for acoustics. There is only
one parameter in this problem: the forcing frequency at the horn inlet.
We impose a radiation boundary condition on the outer boundary.


<br><br>In the Online stage we write out the reflection coefficient as a function
of frequency. The reflection coefficient gives the ratio of the reflected
wave to the applied wave at the inlet.


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
        #include "libmesh/libmesh.h"
        #include "libmesh/mesh.h"
        #include "libmesh/mesh_generation.h"
        #include "libmesh/gmv_io.h"
        #include "libmesh/equation_systems.h"
        #include "libmesh/dof_map.h"
        #include "libmesh/getpot.h"
        #include "libmesh/elem.h"
        
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
        
        #ifndef LIBMESH_USE_COMPLEX_NUMBERS
          libmesh_example_assert(false, "--enable-complex");
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
Parse the input file (reduced_basis_ex7.in) using GetPot
</div>

<div class ="fragment">
<pre>
          std::string parameters_filename = "reduced_basis_ex7.in";
          GetPot infile(parameters_filename);
        
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
          mesh.read("horn.msh");
        
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
            equation_systems.add_system&lt;SimpleRBConstruction&gt; ("Acoustics");
        
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
Read in online_N and initialize online parameters
</div>

<div class ="fragment">
<pre>
            Real online_frequency = infile("online_frequency", 0.);
            RBParameters online_mu;
            online_mu.set_value("frequency", online_frequency);
            rb_eval.set_parameters(online_mu);
            rb_eval.print_parameters();
        
</pre>
</div>
<div class = "comment">
Now do the Online solve using the precomputed reduced basis
</div>

<div class ="fragment">
<pre>
            rb_eval.rb_solve(rb_eval.get_n_basis_functions());
        
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
              
              GMVIO(mesh).write_equation_systems ("RB_sol.gmv",equation_systems);
            }
            
</pre>
</div>
<div class = "comment">
Now do a sweep over frequencies and write out the reflection coefficient for each frequency
</div>

<div class ="fragment">
<pre>
            std::ofstream reflection_coeffs_out("reflection_coefficients.dat");
            
            Real n_frequencies = infile("n_frequencies", 0.);
            Real delta_f = (rb_eval.get_parameter_max("frequency") - rb_eval.get_parameter_min("frequency")) / (n_frequencies-1);
            for(unsigned int freq_i=0; freq_i&lt;n_frequencies; freq_i++)
            {
              Real frequency = rb_eval.get_parameter_min("frequency") + freq_i * delta_f;
              online_mu.set_value("frequency", frequency);
              rb_eval.set_parameters(online_mu);
              rb_eval.rb_solve(rb_eval.get_n_basis_functions());
              
              Number complex_one(1., 0.);
              reflection_coeffs_out &lt;&lt; frequency &lt;&lt; " " &lt;&lt; std::abs(rb_eval.RB_outputs[0] - complex_one) &lt;&lt; std::endl;
            }
            reflection_coeffs_out.close();
        
          }
        
        #endif
        
          return 0;
        }
        
</pre>
</div>

<a name="nocomments"></a> 
<br><br><br> <h1> The source file assembly.h without comments: </h1> 
<pre> 
  #ifndef __assembly_h__
  #define __assembly_h__
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/sparse_matrix.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/numeric_vector.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/dense_matrix.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/dense_vector.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/fe.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/fe_interface.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/fe_base.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/elem_assembly.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/quadrature_gauss.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/rb_theta.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/rb_assembly_expansion.h&quot;</FONT></B>
  
  #define damping_epsilon 0.001
  #define R_rad 12.0
  
  #ifdef LIBMESH_USE_COMPLEX_NUMBERS
  
  using libMesh::ElemAssembly;
  using libMesh::FEInterface;
  using libMesh::FEMContext;
  using libMesh::Number;
  using libMesh::Point;
  using libMesh::RBAssemblyExpansion;
  using libMesh::RBParameters;
  using libMesh::RBTheta;
  using libMesh::RBThetaExpansion;
  using libMesh::Real;
  using libMesh::RealGradient;
  using libMesh::MeshBase;
  using libMesh::libmesh_conj;
  
  <B><FONT COLOR="#228B22">struct</FONT></B> ThetaA0 : RBTheta { <B><FONT COLOR="#228B22">virtual</FONT></B> Number evaluate(<B><FONT COLOR="#228B22">const</FONT></B> RBParameters&amp; mu) { Number val(1., mu.get_value(<B><FONT COLOR="#BC8F8F">&quot;frequency&quot;</FONT></B>)*damping_epsilon); <B><FONT COLOR="#A020F0">return</FONT></B> val; } };
  <B><FONT COLOR="#228B22">struct</FONT></B> ThetaA1 : RBTheta { <B><FONT COLOR="#228B22">virtual</FONT></B> Number evaluate(<B><FONT COLOR="#228B22">const</FONT></B> RBParameters&amp; mu) { Number val(-mu.get_value(<B><FONT COLOR="#BC8F8F">&quot;frequency&quot;</FONT></B>)*mu.get_value(<B><FONT COLOR="#BC8F8F">&quot;frequency&quot;</FONT></B>), 0.); <B><FONT COLOR="#A020F0">return</FONT></B> val; } };
  <B><FONT COLOR="#228B22">struct</FONT></B> ThetaA2 : RBTheta { <B><FONT COLOR="#228B22">virtual</FONT></B> Number evaluate(<B><FONT COLOR="#228B22">const</FONT></B> RBParameters&amp; mu) { Number val(0., mu.get_value(<B><FONT COLOR="#BC8F8F">&quot;frequency&quot;</FONT></B>)); <B><FONT COLOR="#A020F0">return</FONT></B> val; } };
  <B><FONT COLOR="#228B22">struct</FONT></B> ThetaA3 : RBTheta { <B><FONT COLOR="#228B22">virtual</FONT></B> Number evaluate(<B><FONT COLOR="#228B22">const</FONT></B> RBParameters&amp; mu) { Number val(0.5/R_rad, mu.get_value(<B><FONT COLOR="#BC8F8F">&quot;frequency&quot;</FONT></B>)); <B><FONT COLOR="#A020F0">return</FONT></B> val; } };
  
  <B><FONT COLOR="#228B22">struct</FONT></B> ThetaF0 : RBTheta { <B><FONT COLOR="#228B22">virtual</FONT></B> Number evaluate(<B><FONT COLOR="#228B22">const</FONT></B> RBParameters&amp; mu) { Number val(0., 2.*mu.get_value(<B><FONT COLOR="#BC8F8F">&quot;frequency&quot;</FONT></B>)); <B><FONT COLOR="#A020F0">return</FONT></B> val; } };
  
  <B><FONT COLOR="#228B22">struct</FONT></B> ThetaOutput0 : RBTheta { <B><FONT COLOR="#228B22">virtual</FONT></B> Number evaluate(<B><FONT COLOR="#228B22">const</FONT></B> RBParameters&amp; ) { Number val(1., 0.); <B><FONT COLOR="#A020F0">return</FONT></B> val; } };
  
  <B><FONT COLOR="#228B22">struct</FONT></B> AcousticsInnerProduct : ElemAssembly
  {
    <B><FONT COLOR="#228B22">virtual</FONT></B> <B><FONT COLOR="#228B22">void</FONT></B> interior_assembly(FEMContext &amp;c)
    {
      <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> p_var = 0;
  
      <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Real&gt; &amp;JxW =
        c.element_fe_var[p_var]-&gt;get_JxW();
  
      <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;Real&gt; &gt;&amp; phi =
        c.element_fe_var[p_var]-&gt;get_phi();
        
      <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;RealGradient&gt; &gt;&amp; dphi =
        c.element_fe_var[p_var]-&gt;get_dphi();
  
      <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_p_dofs = c.dof_indices_var[p_var].size();
  
      <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_qpoints = (c.get_element_qrule())-&gt;n_points();
  
      <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> qp=0; qp != n_qpoints; qp++)
        <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i != n_p_dofs; i++)
          <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j=0; j != n_p_dofs; j++)
            c.elem_jacobian(i,j) += JxW[qp] * (dphi[j][qp](0)*libmesh_conj(dphi[i][qp](0)) +
                                               dphi[j][qp](1)*libmesh_conj(dphi[i][qp](1)) + 
                                               (phi[j][qp]*libmesh_conj(phi[i][qp])) );
    }
  };
  
  <B><FONT COLOR="#228B22">struct</FONT></B> A0 : ElemAssembly
  {
    <B><FONT COLOR="#228B22">virtual</FONT></B> <B><FONT COLOR="#228B22">void</FONT></B> interior_assembly(FEMContext &amp;c)
    {
      <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> p_var = 0;
  
      <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Real&gt; &amp;JxW =
        c.element_fe_var[p_var]-&gt;get_JxW();
  
      <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;RealGradient&gt; &gt;&amp; dphi =
        c.element_fe_var[p_var]-&gt;get_dphi();
  
      <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_p_dofs = c.dof_indices_var[p_var].size();
  
      <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_qpoints = (c.get_element_qrule())-&gt;n_points();
  
      <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> qp=0; qp != n_qpoints; qp++)
        <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i != n_p_dofs; i++)
          <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j=0; j != n_p_dofs; j++)
            c.elem_jacobian(i,j) += JxW[qp] * (dphi[j][qp](0)*libmesh_conj(dphi[i][qp](0)) +
                                               dphi[j][qp](1)*libmesh_conj(dphi[i][qp](1)));
    }
  };
  
  
  <B><FONT COLOR="#228B22">struct</FONT></B> A1 : ElemAssembly
  {
    <B><FONT COLOR="#228B22">virtual</FONT></B> <B><FONT COLOR="#228B22">void</FONT></B> interior_assembly(FEMContext &amp;c)
    {
      <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> p_var = 0;
  
      <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Real&gt; &amp;JxW =
        c.element_fe_var[p_var]-&gt;get_JxW();
  
      <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;Real&gt; &gt;&amp; phi =
        c.element_fe_var[p_var]-&gt;get_phi();
  
      <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_p_dofs = c.dof_indices_var[p_var].size();
  
      <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_qpoints = (c.get_element_qrule())-&gt;n_points();
  
      <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> qp=0; qp != n_qpoints; qp++)
        <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i != n_p_dofs; i++)
          <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j=0; j != n_p_dofs; j++)
            c.elem_jacobian(i,j) += JxW[qp] * (phi[j][qp]*libmesh_conj(phi[i][qp]));
    }
  };
  
  <B><FONT COLOR="#228B22">struct</FONT></B> A2 : ElemAssembly
  {
    <B><FONT COLOR="#228B22">virtual</FONT></B> <B><FONT COLOR="#228B22">void</FONT></B> boundary_assembly(FEMContext &amp;c)
    {
      <B><FONT COLOR="#A020F0">if</FONT></B>( mesh-&gt;boundary_info-&gt;has_boundary_id (c.elem, c.side, 1) ) <I><FONT COLOR="#B22222">// Forcing on the horn &quot;inlet&quot;
</FONT></I>      {
        <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> p_var = 0;
  
        <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Real&gt; &amp;JxW_face =
          c.side_fe_var[p_var]-&gt;get_JxW();
  
        <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;Real&gt; &gt;&amp; phi_face =
          c.side_fe_var[p_var]-&gt;get_phi();
  
        <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_p_dofs = c.dof_indices_var[p_var].size();
  
        <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_sidepoints = (c.side_qrule)-&gt;n_points();
  
        <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> qp=0; qp != n_sidepoints; qp++)
          <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i != n_p_dofs; i++)
            <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j=0; j != n_p_dofs; j++)
              c.elem_jacobian(i,j) += JxW_face[qp] * phi_face[j][qp] * libmesh_conj(phi_face[i][qp]);
      }
    }
    
    MeshBase* mesh;
  };
  
  <B><FONT COLOR="#228B22">struct</FONT></B> A3 : ElemAssembly
  {
    <B><FONT COLOR="#228B22">virtual</FONT></B> <B><FONT COLOR="#228B22">void</FONT></B> boundary_assembly(FEMContext &amp;c)
    {
      <B><FONT COLOR="#A020F0">if</FONT></B>( mesh-&gt;boundary_info-&gt;has_boundary_id (c.elem, c.side, 2) ) <I><FONT COLOR="#B22222">// Radiation condition on the &quot;bubble&quot;
</FONT></I>      {
        <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> p_var = 0;
  
        <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Real&gt; &amp;JxW_face =
          c.side_fe_var[p_var]-&gt;get_JxW();
  
        <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;Real&gt; &gt;&amp; phi_face =
          c.side_fe_var[p_var]-&gt;get_phi();
  
        <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_p_dofs = c.dof_indices_var[p_var].size();
  
        <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_sidepoints = (c.side_qrule)-&gt;n_points();
  
        <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> qp=0; qp != n_sidepoints; qp++)
          <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i != n_p_dofs; i++)
            <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j=0; j != n_p_dofs; j++)
              c.elem_jacobian(i,j) += JxW_face[qp] * phi_face[j][qp] * libmesh_conj(phi_face[i][qp]);
      }
    }
    
    MeshBase* mesh;
  };
  
  <B><FONT COLOR="#228B22">struct</FONT></B> F0 : ElemAssembly
  {
    <B><FONT COLOR="#228B22">virtual</FONT></B> <B><FONT COLOR="#228B22">void</FONT></B> boundary_assembly(FEMContext &amp;c)
    {
      <B><FONT COLOR="#A020F0">if</FONT></B>( mesh-&gt;boundary_info-&gt;has_boundary_id (c.elem, c.side, 1) ) <I><FONT COLOR="#B22222">// Output is calculated on the horn &quot;inlet&quot;
</FONT></I>      {
        <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> p_var = 0;
  
        <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Real&gt; &amp;JxW_face =
          c.side_fe_var[p_var]-&gt;get_JxW();
  
        <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;Real&gt; &gt;&amp; phi_face =
          c.side_fe_var[p_var]-&gt;get_phi();
  
        <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_p_dofs = c.dof_indices_var[p_var].size();
  
        <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_sidepoints = (c.side_qrule)-&gt;n_points();
  
        <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> qp=0; qp != n_sidepoints; qp++)
          <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i != n_p_dofs; i++)
              c.elem_residual(i) += JxW_face[qp] * libmesh_conj(phi_face[i][qp]);
      }
    }
    
    MeshBase* mesh;
  };
  
  <B><FONT COLOR="#228B22">struct</FONT></B> Output0 : ElemAssembly
  {
    <B><FONT COLOR="#228B22">virtual</FONT></B> <B><FONT COLOR="#228B22">void</FONT></B> boundary_assembly(FEMContext &amp;c)
    {
      <B><FONT COLOR="#A020F0">if</FONT></B>( mesh-&gt;boundary_info-&gt;has_boundary_id (c.elem, c.side, 1) ) <I><FONT COLOR="#B22222">// Forcing on the horn &quot;inlet&quot;
</FONT></I>      {
        <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> p_var = 0;
  
        <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Real&gt; &amp;JxW_face =
          c.side_fe_var[p_var]-&gt;get_JxW();
  
        <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;Real&gt; &gt;&amp; phi_face =
          c.side_fe_var[p_var]-&gt;get_phi();
  
        <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_p_dofs = c.dof_indices_var[p_var].size();
  
        <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_sidepoints = (c.side_qrule)-&gt;n_points();
  
        <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> qp=0; qp != n_sidepoints; qp++)
          <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i != n_p_dofs; i++)
              c.elem_residual(i) += JxW_face[qp] * libmesh_conj(phi_face[i][qp]);
      }
    }
    
    MeshBase* mesh;
  };
  
  <B><FONT COLOR="#228B22">struct</FONT></B> AcousticsRBThetaExpansion : RBThetaExpansion
  {
  
    <I><FONT COLOR="#B22222">/**
     * Constructor.
     */</FONT></I>
    AcousticsRBThetaExpansion()
    {
      attach_A_theta(&amp;theta_a_0);   <I><FONT COLOR="#B22222">// Attach the lhs theta
</FONT></I>      attach_A_theta(&amp;theta_a_1);
      attach_A_theta(&amp;theta_a_2);
      attach_A_theta(&amp;theta_a_3);
  
      attach_F_theta(&amp;theta_f_0);    <I><FONT COLOR="#B22222">// Attach the rhs theta
</FONT></I>      
      attach_output_theta(&amp;theta_output_0);
    }
  
    ThetaA0 theta_a_0;
    ThetaA1 theta_a_1;
    ThetaA2 theta_a_2;
    ThetaA3 theta_a_3;
    ThetaF0 theta_f_0;
    ThetaOutput0 theta_output_0;
  };
  
  <B><FONT COLOR="#228B22">struct</FONT></B> AcousticsRBAssemblyExpansion : RBAssemblyExpansion
  {
  
    <I><FONT COLOR="#B22222">/**
     * Constructor.
     */</FONT></I>
    AcousticsRBAssemblyExpansion(MeshBase&amp; mesh_in)
    {
      A2_assembly.mesh = &amp;mesh_in;
      A3_assembly.mesh = &amp;mesh_in;
      F0_assembly.mesh = &amp;mesh_in;
      Output0_assembly.mesh = &amp;mesh_in;
      
      attach_A_assembly(&amp;A0_assembly); <I><FONT COLOR="#B22222">// Attach the lhs assembly
</FONT></I>      attach_A_assembly(&amp;A1_assembly);
      attach_A_assembly(&amp;A2_assembly);
      attach_A_assembly(&amp;A3_assembly);
      
      attach_F_assembly(&amp;F0_assembly); <I><FONT COLOR="#B22222">// Attach the rhs assembly
</FONT></I>      
      attach_output_assembly(&amp;Output0_assembly);
    }
  
    A0 A0_assembly;
    A1 A1_assembly;
    A2 A2_assembly;
    A3 A3_assembly;
    F0 F0_assembly;
    Output0 Output0_assembly;
    AcousticsInnerProduct acoustics_inner_product;
  };
  
  #endif
  
  #endif
</pre> 
<a name="nocomments"></a> 
<br><br><br> <h1> The source file rb_classes.h without comments: </h1> 
<pre> 
    
    
  
  #ifndef __rb_classes_h__
  #define __rb_classes_h__
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/rb_construction.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/fe_base.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;assembly.h&quot;</FONT></B>
  
  #ifdef LIBMESH_USE_COMPLEX_NUMBERS
  
  using libMesh::AutoPtr;
  using libMesh::DirichletBoundary;
  using libMesh::EquationSystems;
  using libMesh::FEMContext;
  using libMesh::RBConstruction;
  using libMesh::RBEvaluation;
  using libMesh::Real;
  
  
  <B><FONT COLOR="#228B22">class</FONT></B> SimpleRBEvaluation : <B><FONT COLOR="#228B22">public</FONT></B> RBEvaluation
  {
  <B><FONT COLOR="#228B22">public</FONT></B>:
  
    <I><FONT COLOR="#B22222">/**
     * Constructor. Just set the theta expansion.
     */</FONT></I>
    SimpleRBEvaluation()
    {
      set_rb_theta_expansion(acoustics_rb_theta_expansion);
    }
  
    <I><FONT COLOR="#B22222">/**
     * Return a &quot;dummy&quot; stability lower bound factor (for a rigorous error bound
     * this should be a lower bound for the frequency dependent inf-sup constant)
     */</FONT></I>
    <B><FONT COLOR="#228B22">virtual</FONT></B> Real get_stability_lower_bound() { <B><FONT COLOR="#A020F0">return</FONT></B> 1.; }
  
    <I><FONT COLOR="#B22222">/**
     * The object that stores the &quot;theta&quot; expansion of the parameter dependent PDE,
     * i.e. the set of parameter-dependent functions in the affine expansion of the PDE.
     */</FONT></I>
    AcousticsRBThetaExpansion acoustics_rb_theta_expansion;
  
  };
  
  <B><FONT COLOR="#228B22">class</FONT></B> SimpleRBConstruction : <B><FONT COLOR="#228B22">public</FONT></B> RBConstruction
  {
  <B><FONT COLOR="#228B22">public</FONT></B>:
  
    SimpleRBConstruction (EquationSystems&amp; es,
                          <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp; name,
                          <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> number)
    : Parent(es, name, number)
    {}
  
    <I><FONT COLOR="#B22222">/**
     * Destructor.
     */</FONT></I>
    <B><FONT COLOR="#228B22">virtual</FONT></B> ~SimpleRBConstruction () { <B><FONT COLOR="#A020F0">delete</FONT></B> acoustics_rb_assembly_expansion; }
  
    <I><FONT COLOR="#B22222">/**
     * The type of system.
     */</FONT></I>
    <B><FONT COLOR="#228B22">typedef</FONT></B> SimpleRBConstruction sys_type;
  
    <I><FONT COLOR="#B22222">/**
     * The type of the parent.
     */</FONT></I>
    <B><FONT COLOR="#228B22">typedef</FONT></B> RBConstruction Parent;
  
    <I><FONT COLOR="#B22222">/**
     * Initialize data structures.
     */</FONT></I>
    <B><FONT COLOR="#228B22">virtual</FONT></B> <B><FONT COLOR="#228B22">void</FONT></B> init_data()
    {
      p_var = <B><FONT COLOR="#A020F0">this</FONT></B>-&gt;add_variable (<B><FONT COLOR="#BC8F8F">&quot;p&quot;</FONT></B>, SECOND);
  
      <B><FONT COLOR="#5F9EA0">Parent</FONT></B>::init_data();
      
      acoustics_rb_assembly_expansion = <B><FONT COLOR="#A020F0">new</FONT></B> AcousticsRBAssemblyExpansion(get_mesh());
  
      set_rb_assembly_expansion(*acoustics_rb_assembly_expansion);
  
      set_inner_product_assembly(acoustics_rb_assembly_expansion-&gt;acoustics_inner_product);
    }
  
    <I><FONT COLOR="#B22222">/**
     * Pre-request all relevant element data.
     */</FONT></I>
    <B><FONT COLOR="#228B22">virtual</FONT></B> <B><FONT COLOR="#228B22">void</FONT></B> init_context(FEMContext &amp;c)
    {
      c.element_fe_var[p_var]-&gt;get_JxW();
      c.element_fe_var[p_var]-&gt;get_phi();
      c.element_fe_var[p_var]-&gt;get_dphi();
    }
  
    <I><FONT COLOR="#B22222">/**
     * Variable number for pd.
     */</FONT></I>
    <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> p_var;
    
    <I><FONT COLOR="#B22222">/**
     * The object that stores the &quot;assembly&quot; expansion of the parameter dependent PDE,
     * i.e. the objects that define how to assemble the set of parameter-independent
     * operators in the affine expansion of the PDE.
     */</FONT></I>
    AcousticsRBAssemblyExpansion* acoustics_rb_assembly_expansion;
  
  };
  
  #endif
  
  #endif
</pre> 
<a name="nocomments"></a> 
<br><br><br> <h1> The source file reduced_basis_ex7.C without comments: </h1> 
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
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/libmesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/mesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/mesh_generation.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/gmv_io.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/equation_systems.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/dof_map.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/getpot.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/elem.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;rb_classes.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;assembly.h&quot;</FONT></B>
  
  using namespace libMesh;
  
  <B><FONT COLOR="#228B22">int</FONT></B> main (<B><FONT COLOR="#228B22">int</FONT></B> argc, <B><FONT COLOR="#228B22">char</FONT></B>** argv)
  {
    LibMeshInit init (argc, argv);
  
  #ifndef LIBMESH_USE_COMPLEX_NUMBERS
    libmesh_example_assert(false, <B><FONT COLOR="#BC8F8F">&quot;--enable-complex&quot;</FONT></B>);
  #<B><FONT COLOR="#A020F0">else</FONT></B>
  
  #<B><FONT COLOR="#A020F0">if</FONT></B> !defined(LIBMESH_HAVE_XDR)
    libmesh_example_assert(false, <B><FONT COLOR="#BC8F8F">&quot;--enable-xdr&quot;</FONT></B>);
  #elif defined(LIBMESH_DEFAULT_SINGLE_PRECISION)
    libmesh_example_assert(false, <B><FONT COLOR="#BC8F8F">&quot;--disable-singleprecision&quot;</FONT></B>);
  #endif
    libmesh_example_assert(libMesh::default_solver_package() == PETSC_SOLVERS, <B><FONT COLOR="#BC8F8F">&quot;--enable-petsc&quot;</FONT></B>);
  
    libmesh_example_assert(2 &lt;= LIBMESH_DIM, <B><FONT COLOR="#BC8F8F">&quot;2D support&quot;</FONT></B>);
    
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::string parameters_filename = <B><FONT COLOR="#BC8F8F">&quot;reduced_basis_ex7.in&quot;</FONT></B>;
    GetPot infile(parameters_filename);
  
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> dim = 2;                      <I><FONT COLOR="#B22222">// The number of spatial dimensions
</FONT></I>    
    <B><FONT COLOR="#228B22">bool</FONT></B> store_basis_functions = infile(<B><FONT COLOR="#BC8F8F">&quot;store_basis_functions&quot;</FONT></B>, true); <I><FONT COLOR="#B22222">// Do we write the RB basis functions to disk?
</FONT></I>  
    GetPot command_line (argc, argv);
    <B><FONT COLOR="#228B22">int</FONT></B> online_mode = 0;
    <B><FONT COLOR="#A020F0">if</FONT></B> ( command_line.search(1, <B><FONT COLOR="#BC8F8F">&quot;-online_mode&quot;</FONT></B>) )
      online_mode = command_line.next(online_mode);
  
    Mesh mesh (dim);
    mesh.read(<B><FONT COLOR="#BC8F8F">&quot;horn.msh&quot;</FONT></B>);
  
    EquationSystems equation_systems (mesh);
  
    SimpleRBConstruction &amp; rb_con =
      equation_systems.add_system&lt;SimpleRBConstruction&gt; (<B><FONT COLOR="#BC8F8F">&quot;Acoustics&quot;</FONT></B>);
  
    equation_systems.init ();
  
    equation_systems.print_info();
    mesh.print_info();
  
    SimpleRBEvaluation rb_eval;
    
    rb_con.set_rb_evaluation(rb_eval);
  
    <B><FONT COLOR="#A020F0">if</FONT></B>(!online_mode) <I><FONT COLOR="#B22222">// Perform the Offline stage of the RB method
</FONT></I>    {
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
</FONT></I>    {
      rb_eval.read_offline_data_from_files();
      
      Real online_frequency = infile(<B><FONT COLOR="#BC8F8F">&quot;online_frequency&quot;</FONT></B>, 0.);
      RBParameters online_mu;
      online_mu.set_value(<B><FONT COLOR="#BC8F8F">&quot;frequency&quot;</FONT></B>, online_frequency);
      rb_eval.set_parameters(online_mu);
      rb_eval.print_parameters();
  
      rb_eval.rb_solve(rb_eval.get_n_basis_functions());
  
      <B><FONT COLOR="#A020F0">if</FONT></B>(store_basis_functions)
      {
        rb_eval.read_in_basis_functions(rb_con);
        
        rb_con.load_rb_solution();
        
        GMVIO(mesh).write_equation_systems (<B><FONT COLOR="#BC8F8F">&quot;RB_sol.gmv&quot;</FONT></B>,equation_systems);
      }
      
      <B><FONT COLOR="#5F9EA0">std</FONT></B>::ofstream reflection_coeffs_out(<B><FONT COLOR="#BC8F8F">&quot;reflection_coefficients.dat&quot;</FONT></B>);
      
      Real n_frequencies = infile(<B><FONT COLOR="#BC8F8F">&quot;n_frequencies&quot;</FONT></B>, 0.);
      Real delta_f = (rb_eval.get_parameter_max(<B><FONT COLOR="#BC8F8F">&quot;frequency&quot;</FONT></B>) - rb_eval.get_parameter_min(<B><FONT COLOR="#BC8F8F">&quot;frequency&quot;</FONT></B>)) / (n_frequencies-1);
      <B><FONT COLOR="#A020F0">for</FONT></B>(<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> freq_i=0; freq_i&lt;n_frequencies; freq_i++)
      {
        Real frequency = rb_eval.get_parameter_min(<B><FONT COLOR="#BC8F8F">&quot;frequency&quot;</FONT></B>) + freq_i * delta_f;
        online_mu.set_value(<B><FONT COLOR="#BC8F8F">&quot;frequency&quot;</FONT></B>, frequency);
        rb_eval.set_parameters(online_mu);
        rb_eval.rb_solve(rb_eval.get_n_basis_functions());
        
        Number complex_one(1., 0.);
        reflection_coeffs_out &lt;&lt; frequency &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; &quot;</FONT></B> &lt;&lt; std::abs(rb_eval.RB_outputs[0] - complex_one) &lt;&lt; std::endl;
      }
      reflection_coeffs_out.close();
  
    }
  
  #endif
  
    <B><FONT COLOR="#A020F0">return</FONT></B> 0;
  }
  
</pre> 
<a name="output"></a> 
<br><br><br> <h1> The console output of the program: </h1> 
<pre>
***************************************************************
* Running Example reduced_basis_ex7:
*  mpirun -np 12 example-devel -online_mode 0 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
***************************************************************
 
Assertion `false' failed.  Configuring libMesh with --enable-complex may be required to run this code.
************************************************************************************************************************
***             WIDEN YOUR WINDOW TO 120 CHARACTERS.  Use 'enscript -r -fCourier9' to print this document            ***
************************************************************************************************************************

---------------------------------------------- PETSc Performance Summary: ----------------------------------------------

/workspace/libmesh/examples/reduced_basis/reduced_basis_ex7/.libs/lt-example-devel on a intel-12. named hbar.ices.utexas.edu with 12 processors, by benkirk Thu Jan 31 22:20:35 2013
Using Petsc Release Version 3.3.0, Patch 2, Fri Jul 13 15:42:00 CDT 2012 

                         Max       Max/Min        Avg      Total 
Time (sec):           2.829e-03      1.00101   2.828e-03
Objects:              1.000e+00      1.00000   1.000e+00
Flops:                0.000e+00      0.00000   0.000e+00  0.000e+00
Flops/sec:            0.000e+00      0.00000   0.000e+00  0.000e+00
MPI Messages:         0.000e+00      0.00000   0.000e+00  0.000e+00
MPI Message Lengths:  0.000e+00      0.00000   0.000e+00  0.000e+00
MPI Reductions:       1.000e+00      1.00000

Flop counting convention: 1 flop = 1 real number operation of type (multiply/divide/add/subtract)
                            e.g., VecAXPY() for real vectors of length N --> 2N flops
                            and VecAXPY() for complex vectors of length N --> 8N flops

Summary of Stages:   ----- Time ------  ----- Flops -----  --- Messages ---  -- Message Lengths --  -- Reductions --
                        Avg     %Total     Avg     %Total   counts   %Total     Avg         %Total   counts   %Total 
 0:      Main Stage: 2.7758e-03  98.1%  0.0000e+00   0.0%  0.000e+00   0.0%  0.000e+00        0.0%  0.000e+00   0.0% 

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
      %T - percent time in this phase         %f - percent flops in this phase
      %M - percent messages in this phase     %L - percent message lengths in this phase
      %R - percent reductions in this phase
   Total Mflop/s: 10e-6 * (sum of flops over all processors)/(max time over all processors)
------------------------------------------------------------------------------------------------------------------------
Event                Count      Time (sec)     Flops                             --- Global ---  --- Stage ---   Total
                   Max Ratio  Max     Ratio   Max  Ratio  Mess   Avg len Reduct  %T %f %M %L %R  %T %f %M %L %R Mflop/s
------------------------------------------------------------------------------------------------------------------------

--- Event Stage 0: Main Stage

------------------------------------------------------------------------------------------------------------------------

Memory usage is given in bytes:

Object Type          Creations   Destructions     Memory  Descendants' Mem.
Reports information only for process 0.

--- Event Stage 0: Main Stage

              Viewer     1              0            0     0
========================================================================================================================
Average time to get PetscTime(): 0
Average time for MPI_Barrier(): 1.06335e-05
Average time for zero size MPI_Send(): 1.659e-05
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
sizeof(short) 2 sizeof(int) 4 sizeof(long) 8 sizeof(void*) 8 sizeof(PetscScalar) 8 sizeof(PetscInt) 4
Configure run at: Thu Nov  8 11:21:02 2012
Configure options: --with-debugging=false --COPTFLAGS=-O3 --CXXOPTFLAGS=-O3 --FOPTFLAGS=-O3 --with-clanguage=C++ --with-shared-libraries=1 --with-mpi-dir=/opt/apps/ossw/libraries/mpich2/mpich2-1.4.1p1/sl6/intel-12.1 --with-mumps=true --download-mumps=1 --with-metis=true --download-metis=1 --with-parmetis=true --download-parmetis=1 --with-superlu=true --download-superlu=1 --with-superludir=true --download-superlu_dist=1 --with-blacs=true --download-blacs=1 --with-scalapack=true --download-scalapack=1 --with-hypre=true --download-hypre=1 --with-blas-lib="[/opt/apps/sysnet/intel/12.1/mkl/10.3.12.361/lib/intel64/libmkl_intel_lp64.so,/opt/apps/sysnet/intel/12.1/mkl/10.3.12.361/lib/intel64/libmkl_sequential.so,/opt/apps/sysnet/intel/12.1/mkl/10.3.12.361/lib/intel64/libmkl_core.so]" --with-lapack-lib="[/opt/apps/sysnet/intel/12.1/mkl/10.3.12.361/lib/intel64/libmkl_lapack95_lp64.a]"
-----------------------------------------
Libraries compiled on Thu Nov  8 11:21:02 2012 on daedalus.ices.utexas.edu 
Machine characteristics: Linux-2.6.32-279.1.1.el6.x86_64-x86_64-with-redhat-6.3-Carbon
Using PETSc directory: /opt/apps/ossw/libraries/petsc/petsc-3.3-p2
Using PETSc arch: intel-12.1-mkl-intel-10.3.12.361-mpich2-1.4.1p1-cxx-opt
-----------------------------------------

Using C compiler: /opt/apps/ossw/libraries/mpich2/mpich2-1.4.1p1/sl6/intel-12.1/bin/mpicxx  -wd1572 -O3   -fPIC   ${COPTFLAGS} ${CFLAGS}
Using Fortran compiler: /opt/apps/ossw/libraries/mpich2/mpich2-1.4.1p1/sl6/intel-12.1/bin/mpif90  -fPIC -O3   ${FOPTFLAGS} ${FFLAGS} 
-----------------------------------------

Using include paths: -I/opt/apps/ossw/libraries/petsc/petsc-3.3-p2/intel-12.1-mkl-intel-10.3.12.361-mpich2-1.4.1p1-cxx-opt/include -I/opt/apps/ossw/libraries/petsc/petsc-3.3-p2/include -I/opt/apps/ossw/libraries/petsc/petsc-3.3-p2/include -I/opt/apps/ossw/libraries/petsc/petsc-3.3-p2/intel-12.1-mkl-intel-10.3.12.361-mpich2-1.4.1p1-cxx-opt/include -I/opt/apps/ossw/libraries/mpich2/mpich2-1.4.1p1/sl6/intel-12.1/include
-----------------------------------------

Using C linker: /opt/apps/ossw/libraries/mpich2/mpich2-1.4.1p1/sl6/intel-12.1/bin/mpicxx
Using Fortran linker: /opt/apps/ossw/libraries/mpich2/mpich2-1.4.1p1/sl6/intel-12.1/bin/mpif90
Using libraries: -Wl,-rpath,/opt/apps/ossw/libraries/petsc/petsc-3.3-p2/intel-12.1-mkl-intel-10.3.12.361-mpich2-1.4.1p1-cxx-opt/lib -L/opt/apps/ossw/libraries/petsc/petsc-3.3-p2/intel-12.1-mkl-intel-10.3.12.361-mpich2-1.4.1p1-cxx-opt/lib -lpetsc -lX11 -Wl,-rpath,/opt/apps/ossw/libraries/petsc/petsc-3.3-p2/intel-12.1-mkl-intel-10.3.12.361-mpich2-1.4.1p1-cxx-opt/lib -L/opt/apps/ossw/libraries/petsc/petsc-3.3-p2/intel-12.1-mkl-intel-10.3.12.361-mpich2-1.4.1p1-cxx-opt/lib -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord -lHYPRE -lpthread -lsuperlu_dist_3.0 -lparmetis -lmetis -lscalapack -lblacs -lsuperlu_4.3 -Wl,-rpath,/opt/apps/sysnet/intel/12.1/mkl/10.3.12.361/lib/intel64 -L/opt/apps/sysnet/intel/12.1/mkl/10.3.12.361/lib/intel64 -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -Wl,-rpath,/opt/apps/ossw/libraries/mpich2/mpich2-1.4.1p1/sl6/intel-12.1/lib -L/opt/apps/ossw/libraries/mpich2/mpich2-1.4.1p1/sl6/intel-12.1/lib -Wl,-rpath,/opt/apps/sysnet/intel/12.1/composer_xe_2011_sp1.7.256/compiler/lib/intel64 -L/opt/apps/sysnet/intel/12.1/composer_xe_2011_sp1.7.256/compiler/lib/intel64 -Wl,-rpath,/usr/lib/gcc/x86_64-redhat-linux/4.4.6 -L/usr/lib/gcc/x86_64-redhat-linux/4.4.6 -lmpichf90 -lifport -lifcore -lm -lm -lmpichcxx -ldl -lmpich -lopa -lmpl -lrt -lpthread -limf -lsvml -lipgo -ldecimal -lcilkrts -lstdc++ -lgcc_s -lirc -lirc_s -ldl 
-----------------------------------------


 ----------------------------------------------------------------------------------------------------------------------
| Processor id:   0                                                                                                    |
| Num Processors: 12                                                                                                   |
| Time:           Thu Jan 31 22:20:35 2013                                                                             |
| OS:             Linux                                                                                                |
| HostName:       hbar.ices.utexas.edu                                                                                 |
| OS Release:     2.6.32-279.1.1.el6.x86_64                                                                            |
| OS Version:     #1 SMP Tue Jul 10 11:24:23 CDT 2012                                                                  |
| Machine:        x86_64                                                                                               |
| Username:       benkirk                                                                                              |
| Configuration:  ./configure  '--enable-everything'                                                                   |
|  '--prefix=/workspace/libmesh/install'                                                                               |
|  'CXX=icpc'                                                                                                          |
|  'CC=icc'                                                                                                            |
|  'FC=ifort'                                                                                                          |
|  'F77=ifort'                                                                                                         |
|  'PETSC_DIR=/opt/apps/ossw/libraries/petsc/petsc-3.3-p2'                                                             |
|  'PETSC_ARCH=intel-12.1-mkl-intel-10.3.12.361-mpich2-1.4.1p1-cxx-opt'                                                |
|  'SLEPC_DIR=/opt/apps/ossw/libraries/slepc/slepc-3.3-p2-petsc-3.3-p2-cxx-opt'                                        |
|  'TRILINOS_DIR=/opt/apps/ossw/libraries/trilinos/trilinos-10.12.2/sl6/intel-12.1/mpich2-1.4.1p1/mkl-intel-10.3.12.361'|
|  'VTK_DIR=/opt/apps/ossw/libraries/vtk/vtk-5.10.0/sl6/intel-12.1'                                                    |
 ----------------------------------------------------------------------------------------------------------------------
 -----------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=0.02404, Active time=0.000486                                             |
 -----------------------------------------------------------------------------------------------------------
| Event                         nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                         w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|-----------------------------------------------------------------------------------------------------------|
|                                                                                                           |
|                                                                                                           |
| Parallel                                                                                                  |
|   max(scalar)                 8         0.0001      0.000011    0.0001      0.000011    18.52    18.52    |
|   max(vector)                 2         0.0000      0.000022    0.0001      0.000047    8.85     19.34    |
|   min(bool)                   10        0.0001      0.000007    0.0001      0.000007    15.23    15.23    |
|   min(scalar)                 8         0.0002      0.000028    0.0002      0.000028    46.30    46.30    |
|   min(vector)                 2         0.0001      0.000027    0.0001      0.000057    11.11    23.25    |
 -----------------------------------------------------------------------------------------------------------
| Totals:                       30        0.0005                                          100.00            |
 -----------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example reduced_basis_ex7:
*  mpirun -np 12 example-devel -online_mode 0 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
***************************************************************
***************************************************************
* Running Example reduced_basis_ex7:
*  mpirun -np 12 example-devel -online_mode 1 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
***************************************************************
 
Assertion `false' failed.  Configuring libMesh with --enable-complex may be required to run this code.
************************************************************************************************************************
***             WIDEN YOUR WINDOW TO 120 CHARACTERS.  Use 'enscript -r -fCourier9' to print this document            ***
************************************************************************************************************************

---------------------------------------------- PETSc Performance Summary: ----------------------------------------------

/workspace/libmesh/examples/reduced_basis/reduced_basis_ex7/.libs/lt-example-devel on a intel-12. named hbar.ices.utexas.edu with 12 processors, by benkirk Thu Jan 31 22:20:35 2013
Using Petsc Release Version 3.3.0, Patch 2, Fri Jul 13 15:42:00 CDT 2012 

                         Max       Max/Min        Avg      Total 
Time (sec):           4.390e-03      1.50803   3.942e-03
Objects:              1.000e+00      1.00000   1.000e+00
Flops:                0.000e+00      0.00000   0.000e+00  0.000e+00
Flops/sec:            0.000e+00      0.00000   0.000e+00  0.000e+00
MPI Messages:         0.000e+00      0.00000   0.000e+00  0.000e+00
MPI Message Lengths:  0.000e+00      0.00000   0.000e+00  0.000e+00
MPI Reductions:       1.000e+00      1.00000

Flop counting convention: 1 flop = 1 real number operation of type (multiply/divide/add/subtract)
                            e.g., VecAXPY() for real vectors of length N --> 2N flops
                            and VecAXPY() for complex vectors of length N --> 8N flops

Summary of Stages:   ----- Time ------  ----- Flops -----  --- Messages ---  -- Message Lengths --  -- Reductions --
                        Avg     %Total     Avg     %Total   counts   %Total     Avg         %Total   counts   %Total 
 0:      Main Stage: 3.9062e-03  99.1%  0.0000e+00   0.0%  0.000e+00   0.0%  0.000e+00        0.0%  0.000e+00   0.0% 

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
      %T - percent time in this phase         %f - percent flops in this phase
      %M - percent messages in this phase     %L - percent message lengths in this phase
      %R - percent reductions in this phase
   Total Mflop/s: 10e-6 * (sum of flops over all processors)/(max time over all processors)
------------------------------------------------------------------------------------------------------------------------
Event                Count      Time (sec)     Flops                             --- Global ---  --- Stage ---   Total
                   Max Ratio  Max     Ratio   Max  Ratio  Mess   Avg len Reduct  %T %f %M %L %R  %T %f %M %L %R Mflop/s
------------------------------------------------------------------------------------------------------------------------

--- Event Stage 0: Main Stage

------------------------------------------------------------------------------------------------------------------------

Memory usage is given in bytes:

Object Type          Creations   Destructions     Memory  Descendants' Mem.
Reports information only for process 0.

--- Event Stage 0: Main Stage

              Viewer     1              0            0     0
========================================================================================================================
Average time to get PetscTime(): 0
Average time for MPI_Barrier(): 5.89848e-05
Average time for zero size MPI_Send(): 4.18425e-05
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
sizeof(short) 2 sizeof(int) 4 sizeof(long) 8 sizeof(void*) 8 sizeof(PetscScalar) 8 sizeof(PetscInt) 4
Configure run at: Thu Nov  8 11:21:02 2012
Configure options: --with-debugging=false --COPTFLAGS=-O3 --CXXOPTFLAGS=-O3 --FOPTFLAGS=-O3 --with-clanguage=C++ --with-shared-libraries=1 --with-mpi-dir=/opt/apps/ossw/libraries/mpich2/mpich2-1.4.1p1/sl6/intel-12.1 --with-mumps=true --download-mumps=1 --with-metis=true --download-metis=1 --with-parmetis=true --download-parmetis=1 --with-superlu=true --download-superlu=1 --with-superludir=true --download-superlu_dist=1 --with-blacs=true --download-blacs=1 --with-scalapack=true --download-scalapack=1 --with-hypre=true --download-hypre=1 --with-blas-lib="[/opt/apps/sysnet/intel/12.1/mkl/10.3.12.361/lib/intel64/libmkl_intel_lp64.so,/opt/apps/sysnet/intel/12.1/mkl/10.3.12.361/lib/intel64/libmkl_sequential.so,/opt/apps/sysnet/intel/12.1/mkl/10.3.12.361/lib/intel64/libmkl_core.so]" --with-lapack-lib="[/opt/apps/sysnet/intel/12.1/mkl/10.3.12.361/lib/intel64/libmkl_lapack95_lp64.a]"
-----------------------------------------
Libraries compiled on Thu Nov  8 11:21:02 2012 on daedalus.ices.utexas.edu 
Machine characteristics: Linux-2.6.32-279.1.1.el6.x86_64-x86_64-with-redhat-6.3-Carbon
Using PETSc directory: /opt/apps/ossw/libraries/petsc/petsc-3.3-p2
Using PETSc arch: intel-12.1-mkl-intel-10.3.12.361-mpich2-1.4.1p1-cxx-opt
-----------------------------------------

Using C compiler: /opt/apps/ossw/libraries/mpich2/mpich2-1.4.1p1/sl6/intel-12.1/bin/mpicxx  -wd1572 -O3   -fPIC   ${COPTFLAGS} ${CFLAGS}
Using Fortran compiler: /opt/apps/ossw/libraries/mpich2/mpich2-1.4.1p1/sl6/intel-12.1/bin/mpif90  -fPIC -O3   ${FOPTFLAGS} ${FFLAGS} 
-----------------------------------------

Using include paths: -I/opt/apps/ossw/libraries/petsc/petsc-3.3-p2/intel-12.1-mkl-intel-10.3.12.361-mpich2-1.4.1p1-cxx-opt/include -I/opt/apps/ossw/libraries/petsc/petsc-3.3-p2/include -I/opt/apps/ossw/libraries/petsc/petsc-3.3-p2/include -I/opt/apps/ossw/libraries/petsc/petsc-3.3-p2/intel-12.1-mkl-intel-10.3.12.361-mpich2-1.4.1p1-cxx-opt/include -I/opt/apps/ossw/libraries/mpich2/mpich2-1.4.1p1/sl6/intel-12.1/include
-----------------------------------------

Using C linker: /opt/apps/ossw/libraries/mpich2/mpich2-1.4.1p1/sl6/intel-12.1/bin/mpicxx
Using Fortran linker: /opt/apps/ossw/libraries/mpich2/mpich2-1.4.1p1/sl6/intel-12.1/bin/mpif90
Using libraries: -Wl,-rpath,/opt/apps/ossw/libraries/petsc/petsc-3.3-p2/intel-12.1-mkl-intel-10.3.12.361-mpich2-1.4.1p1-cxx-opt/lib -L/opt/apps/ossw/libraries/petsc/petsc-3.3-p2/intel-12.1-mkl-intel-10.3.12.361-mpich2-1.4.1p1-cxx-opt/lib -lpetsc -lX11 -Wl,-rpath,/opt/apps/ossw/libraries/petsc/petsc-3.3-p2/intel-12.1-mkl-intel-10.3.12.361-mpich2-1.4.1p1-cxx-opt/lib -L/opt/apps/ossw/libraries/petsc/petsc-3.3-p2/intel-12.1-mkl-intel-10.3.12.361-mpich2-1.4.1p1-cxx-opt/lib -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord -lHYPRE -lpthread -lsuperlu_dist_3.0 -lparmetis -lmetis -lscalapack -lblacs -lsuperlu_4.3 -Wl,-rpath,/opt/apps/sysnet/intel/12.1/mkl/10.3.12.361/lib/intel64 -L/opt/apps/sysnet/intel/12.1/mkl/10.3.12.361/lib/intel64 -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -Wl,-rpath,/opt/apps/ossw/libraries/mpich2/mpich2-1.4.1p1/sl6/intel-12.1/lib -L/opt/apps/ossw/libraries/mpich2/mpich2-1.4.1p1/sl6/intel-12.1/lib -Wl,-rpath,/opt/apps/sysnet/intel/12.1/composer_xe_2011_sp1.7.256/compiler/lib/intel64 -L/opt/apps/sysnet/intel/12.1/composer_xe_2011_sp1.7.256/compiler/lib/intel64 -Wl,-rpath,/usr/lib/gcc/x86_64-redhat-linux/4.4.6 -L/usr/lib/gcc/x86_64-redhat-linux/4.4.6 -lmpichf90 -lifport -lifcore -lm -lm -lmpichcxx -ldl -lmpich -lopa -lmpl -lrt -lpthread -limf -lsvml -lipgo -ldecimal -lcilkrts -lstdc++ -lgcc_s -lirc -lirc_s -ldl 
-----------------------------------------


 ----------------------------------------------------------------------------------------------------------------------
| Processor id:   0                                                                                                    |
| Num Processors: 12                                                                                                   |
| Time:           Thu Jan 31 22:20:35 2013                                                                             |
| OS:             Linux                                                                                                |
| HostName:       hbar.ices.utexas.edu                                                                                 |
| OS Release:     2.6.32-279.1.1.el6.x86_64                                                                            |
| OS Version:     #1 SMP Tue Jul 10 11:24:23 CDT 2012                                                                  |
| Machine:        x86_64                                                                                               |
| Username:       benkirk                                                                                              |
| Configuration:  ./configure  '--enable-everything'                                                                   |
|  '--prefix=/workspace/libmesh/install'                                                                               |
|  'CXX=icpc'                                                                                                          |
|  'CC=icc'                                                                                                            |
|  'FC=ifort'                                                                                                          |
|  'F77=ifort'                                                                                                         |
|  'PETSC_DIR=/opt/apps/ossw/libraries/petsc/petsc-3.3-p2'                                                             |
|  'PETSC_ARCH=intel-12.1-mkl-intel-10.3.12.361-mpich2-1.4.1p1-cxx-opt'                                                |
|  'SLEPC_DIR=/opt/apps/ossw/libraries/slepc/slepc-3.3-p2-petsc-3.3-p2-cxx-opt'                                        |
|  'TRILINOS_DIR=/opt/apps/ossw/libraries/trilinos/trilinos-10.12.2/sl6/intel-12.1/mpich2-1.4.1p1/mkl-intel-10.3.12.361'|
|  'VTK_DIR=/opt/apps/ossw/libraries/vtk/vtk-5.10.0/sl6/intel-12.1'                                                    |
 ----------------------------------------------------------------------------------------------------------------------
 -----------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=0.02847, Active time=0.002714                                             |
 -----------------------------------------------------------------------------------------------------------
| Event                         nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                         w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|-----------------------------------------------------------------------------------------------------------|
|                                                                                                           |
|                                                                                                           |
| Parallel                                                                                                  |
|   max(scalar)                 8         0.0003      0.000039    0.0003      0.000039    11.39    11.39    |
|   max(vector)                 2         0.0001      0.000042    0.0003      0.000152    3.10     11.20    |
|   min(bool)                   10        0.0004      0.000036    0.0004      0.000036    13.19    13.19    |
|   min(scalar)                 8         0.0018      0.000230    0.0018      0.000230    67.94    67.94    |
|   min(vector)                 2         0.0001      0.000059    0.0004      0.000187    4.38     13.78    |
 -----------------------------------------------------------------------------------------------------------
| Totals:                       30        0.0027                                          100.00            |
 -----------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example reduced_basis_ex7:
*  mpirun -np 12 example-devel -online_mode 1 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
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
