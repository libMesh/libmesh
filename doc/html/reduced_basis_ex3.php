<?php $root=""; ?>
<?php require($root."navigation.php"); ?>
<html>
<head>
  <?php load_style($root); ?>
</head>
 
<body>
 
<?php make_navigation("reduced_basis_ex3",$root)?>
 
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
        #include "libmesh/transient_rb_theta_expansion.h"
        #include "libmesh/transient_rb_assembly_expansion.h"
        
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
        using libMesh::ElemAssembly;
        using libMesh::FEInterface;
        using libMesh::FEMContext;
        using libMesh::Number;
        using libMesh::Point;
        using libMesh::RBParameters;
        using libMesh::RBTheta;
        using libMesh::Real;
        using libMesh::RealGradient;
        using libMesh::TransientRBThetaExpansion;
        using libMesh::TransientRBAssemblyExpansion;
        
</pre>
</div>
<div class = "comment">
Functors for the parameter-dependent part of the affine decomposition of the PDE
The RHS and outputs just require a constant value of 1, so use a default RBTheta object there
</div>

<div class ="fragment">
<pre>
        struct ThetaA0 : RBTheta { virtual Number evaluate(const RBParameters& )   { return 0.05;  } };
        struct ThetaA1 : RBTheta { virtual Number evaluate(const RBParameters& mu) { return mu.get_value("x_vel"); } };
        struct ThetaA2 : RBTheta { virtual Number evaluate(const RBParameters& mu) { return mu.get_value("y_vel"); } };
        
        struct M0 : ElemAssembly
        {
</pre>
</div>
<div class = "comment">
L2 matrix
</div>

<div class ="fragment">
<pre>
          virtual void interior_assembly(FEMContext &c)
          {
            const unsigned int u_var = 0;
        
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
            unsigned int n_qpoints = (c.get_element_qrule())-&gt;n_points();
        
            for (unsigned int qp=0; qp != n_qpoints; qp++)
              for (unsigned int i=0; i != n_u_dofs; i++)
                for (unsigned int j=0; j != n_u_dofs; j++)
                  c.get_elem_jacobian()(i,j) += JxW[qp] *phi[j][qp]*phi[i][qp];
          }
        };
        
        struct A0 : ElemAssembly
        {
</pre>
</div>
<div class = "comment">
Assemble the Laplacian operator
</div>

<div class ="fragment">
<pre>
          virtual void interior_assembly(FEMContext &c)
          {
            const unsigned int u_var = 0;
        
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
            unsigned int n_qpoints = (c.get_element_qrule())-&gt;n_points();
        
            for (unsigned int qp=0; qp != n_qpoints; qp++)
              for (unsigned int i=0; i != n_u_dofs; i++)
                for (unsigned int j=0; j != n_u_dofs; j++)
                  c.get_elem_jacobian()(i,j) += JxW[qp] * dphi[j][qp]*dphi[i][qp];
          }
        };
        
        
        struct A1 : ElemAssembly
        {
</pre>
</div>
<div class = "comment">
Convection in the x-direction
</div>

<div class ="fragment">
<pre>
          virtual void interior_assembly(FEMContext &c)
          {
            const unsigned int u_var = 0;
        
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
            unsigned int n_qpoints = (c.get_element_qrule())-&gt;n_points();
        
            for (unsigned int qp=0; qp != n_qpoints; qp++)
              for (unsigned int i=0; i != n_u_dofs; i++)
                for (unsigned int j=0; j != n_u_dofs; j++)
                  c.get_elem_jacobian()(i,j) += JxW[qp] *dphi[j][qp](0)*phi[i][qp];
          }
        };
        
        struct A2 : ElemAssembly
        {
</pre>
</div>
<div class = "comment">
Convection in the y-direction
</div>

<div class ="fragment">
<pre>
          virtual void interior_assembly(FEMContext &c)
          {
            const unsigned int u_var = 0;
        
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
            unsigned int n_qpoints = (c.get_element_qrule())-&gt;n_points();
        
            for (unsigned int qp=0; qp != n_qpoints; qp++)
              for (unsigned int i=0; i != n_u_dofs; i++)
                for (unsigned int j=0; j != n_u_dofs; j++)
                  c.get_elem_jacobian()(i,j) += JxW[qp] *dphi[j][qp](1)*phi[i][qp];
          }
        };
        
        struct F0 : ElemAssembly
        {
</pre>
</div>
<div class = "comment">
Source term, 1 throughout the domain
</div>

<div class ="fragment">
<pre>
          virtual void interior_assembly(FEMContext &c)
          {
            const unsigned int u_var = 0;
        
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
            unsigned int n_qpoints = (c.get_element_qrule())-&gt;n_points();
        
            for (unsigned int qp=0; qp != n_qpoints; qp++)
              for (unsigned int i=0; i != n_u_dofs; i++)
                c.get_elem_residual()(i) += JxW[qp] * ( 1.*phi[i][qp] );
          }
        };
        
        struct OutputAssembly : ElemAssembly
        {
          OutputAssembly(Real min_x_in, Real max_x_in,
                         Real min_y_in, Real max_y_in)
                        :
                        min_x(min_x_in),
                        max_x(max_x_in),
                        min_y(min_y_in),
                        max_y(max_y_in)
          {}
        
</pre>
</div>
<div class = "comment">
Output: Average value over the region [min_x,max_x]x[min_y,max_y]
</div>

<div class ="fragment">
<pre>
          virtual void interior_assembly(FEMContext &c)
          {
            const unsigned int u_var = 0;
        
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
            unsigned int n_qpoints = (c.get_element_qrule())-&gt;n_points();
            
            Real output_area = (max_x-min_x) * (max_y-min_y);
        
            Point centroid = c.elem-&gt;centroid();
            if( (min_x &lt;= centroid(0)) && (centroid(0) &lt;= max_x) &&
                (min_y &lt;= centroid(1)) && (centroid(1) &lt;= max_y) )
              for (unsigned int qp=0; qp != n_qpoints; qp++)
                for (unsigned int i=0; i != n_u_dofs; i++)
                  c.get_elem_residual()(i) += JxW[qp] * ( 1.*phi[i][qp] ) / output_area;
          }
          
</pre>
</div>
<div class = "comment">
Member variables that define the output region in 2D
</div>

<div class ="fragment">
<pre>
          Real min_x, max_x, min_y, max_y;
        };
        
</pre>
</div>
<div class = "comment">
Define an RBThetaExpansion class for this PDE
</div>

<div class ="fragment">
<pre>
        struct CDRBThetaExpansion : TransientRBThetaExpansion
        {
        
          /**
           * Constructor.
           */
          CDRBThetaExpansion()
          {
</pre>
</div>
<div class = "comment">
set up the RBThetaExpansion object
</div>

<div class ="fragment">
<pre>
            attach_M_theta(&rb_theta);    // Attach the time-derivative theta
        
            attach_A_theta(&theta_a_0);   // Attach the lhs theta
            attach_A_theta(&theta_a_1);
            attach_A_theta(&theta_a_2);
        
            attach_F_theta(&rb_theta);    // Attach the rhs theta
        
            attach_output_theta(&rb_theta); // Attach output 0 theta
            attach_output_theta(&rb_theta); // Attach output 1 theta
            attach_output_theta(&rb_theta); // Attach output 2 theta
            attach_output_theta(&rb_theta); // Attach output 3 theta
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
          RBTheta rb_theta; // Default RBTheta object, just returns 1.
        };
        
</pre>
</div>
<div class = "comment">
Define an RBAssemblyExpansion class for this PDE
</div>

<div class ="fragment">
<pre>
        struct CDRBAssemblyExpansion : TransientRBAssemblyExpansion
        {
        
          /**
           * Constructor.
           */
          CDRBAssemblyExpansion()
            :
            L0(0.72,0.88,0.72,0.88), // We make sure these output regions conform to the mesh
            L1(0.12,0.28,0.72,0.88),
            L2(0.12,0.28,0.12,0.28),
            L3(0.72,0.88,0.12,0.28)
          {
</pre>
</div>
<div class = "comment">
And set up the RBAssemblyExpansion object
</div>

<div class ="fragment">
<pre>
            attach_M_assembly(&M0_assembly); // Attach the time-derivative assembly
            attach_A_assembly(&A0_assembly); // Attach the lhs assembly
            attach_A_assembly(&A1_assembly);
            attach_A_assembly(&A2_assembly);
            
            attach_F_assembly(&F0_assembly); // Attach the rhs assembly
            
            attach_output_assembly(&L0);       // Attach output 0 assembly
            attach_output_assembly(&L1);       // Attach output 1 assembly
            attach_output_assembly(&L2);       // Attach output 2 assembly
            attach_output_assembly(&L3);       // Attach output 3 assembly
          }
        
</pre>
</div>
<div class = "comment">
The ElemAssembly objects
</div>

<div class ="fragment">
<pre>
          M0 M0_assembly;
          A0 A0_assembly;
          A1 A1_assembly;
          A2 A2_assembly;
          F0 F0_assembly;
          OutputAssembly L0;
          OutputAssembly L1;
          OutputAssembly L2;
          OutputAssembly L3;
        };
        
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
        
        #include "libmesh/transient_rb_construction.h"
        #include "libmesh/fe_base.h"
        
</pre>
</div>
<div class = "comment">
local include
</div>

<div class ="fragment">
<pre>
        #include "assembly.h"
        
</pre>
</div>
<div class = "comment">
Bring in bits from the libMesh namespace.
Just the bits we're using, since this is a header.
</div>

<div class ="fragment">
<pre>
        using libMesh::EquationSystems;
        using libMesh::FEMContext;
        using libMesh::RBConstruction;
        using libMesh::RBEvaluation;
        using libMesh::Real;
        using libMesh::TransientRBEvaluation;
        using libMesh::TransientRBConstruction;
        
        
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
        class SimpleRBEvaluation : public TransientRBEvaluation
        {
        public:
        
          /**
           * Constructor. Just set the theta expansion.
           */
          SimpleRBEvaluation()
          {
            set_rb_theta_expansion(cd_rb_theta_expansion);
          }
        
          /**
           * The coercivity constant is bounded below by 0.05.
           */
          virtual Real get_stability_lower_bound() { return 0.05; }
        
          /**
           * The object that stores the "theta" expansion of the parameter dependent PDE,
           * i.e. the set of parameter-dependent functions in the affine expansion of the PDE.
           */
          CDRBThetaExpansion cd_rb_theta_expansion;
        
        };
        
</pre>
</div>
<div class = "comment">
A simple subclass of Construction, which just needs to override build_rb_evaluation
in order to build a SimpleRBEvaluation object, rather than an RBEvaluation object.
</div>

<div class ="fragment">
<pre>
        class SimpleRBConstruction : public TransientRBConstruction
        {
        public:
        
          SimpleRBConstruction (EquationSystems& es,
                                const std::string& name,
                                const unsigned int number)
          : Parent(es, name, number),
            dirichlet_bc(AutoPtr&lt;DirichletBoundary&gt;(NULL))
          {}
        
          /**
           * Destructor.
           */
          virtual ~SimpleRBConstruction () {}
        
          /**
           * The type of system.
           */
          typedef SimpleRBConstruction sys_type;
        
          /**
           * The type of the parent.
           */
          typedef TransientRBConstruction Parent;
        
          /**
           * Initialize data structures.
           */
          virtual void init_data()
          {
            u_var = this-&gt;add_variable ("u", FIRST);
        
</pre>
</div>
<div class = "comment">
Generate a DirichletBoundary object
</div>

<div class ="fragment">
<pre>
            dirichlet_bc = build_zero_dirichlet_boundary_object();
        
</pre>
</div>
<div class = "comment">
Set the Dirichet boundary IDs
and the Dirichlet boundary variable numbers
</div>

<div class ="fragment">
<pre>
            dirichlet_bc-&gt;b.insert(0);
            dirichlet_bc-&gt;b.insert(1);
            dirichlet_bc-&gt;b.insert(2);
            dirichlet_bc-&gt;b.insert(3);
            dirichlet_bc-&gt;variables.push_back(u_var);
            
</pre>
</div>
<div class = "comment">
Attach dirichlet_bc (must do this _before_ Parent::init_data)
</div>

<div class ="fragment">
<pre>
            get_dof_map().add_dirichlet_boundary(*dirichlet_bc);
        
            Parent::init_data();
        
</pre>
</div>
<div class = "comment">
Set the rb_assembly_expansion for this Construction object.
</div>

<div class ="fragment">
<pre>
            set_rb_assembly_expansion(cd_rb_assembly_expansion);
        
</pre>
</div>
<div class = "comment">
We need to define an inner product matrix for this problem
</div>

<div class ="fragment">
<pre>
            set_inner_product_assembly(cd_rb_assembly_expansion.A0_assembly);
        
</pre>
</div>
<div class = "comment">
We also need an "L2 matrix" in order to compute the time-dependent error bound
</div>

<div class ="fragment">
<pre>
            set_L2_assembly(cd_rb_assembly_expansion.M0_assembly);
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
            c.element_fe_var[u_var]-&gt;get_JxW();
            c.element_fe_var[u_var]-&gt;get_phi();
            c.element_fe_var[u_var]-&gt;get_dphi();
          }
        
          /**
           * Variable number for u.
           */
          unsigned int u_var;
          
          /**
           * The object that stores the "assembly" expansion of the parameter dependent PDE,
           * i.e. the objects that define how to assemble the set of parameter-independent
           * operators in the affine expansion of the PDE.
           */
          CDRBAssemblyExpansion cd_rb_assembly_expansion;
        
          /**
           * The object that defines which degrees of freedom are on a Dirichlet boundary.
           */
          AutoPtr&lt;DirichletBoundary&gt; dirichlet_bc;
        
        };
        
        #endif
</pre>
</div>

<a name="comments"></a> 
<br><br><br> <h1> The source file reduced_basis_ex3.C with comments: </h1> 
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
<h1>Reduced Basis Example 3 - Transient Reduced Basis Problem</h1>


<br><br>In this example problem we use the Certified Reduced Basis method
to solve a transient convection-diffusion problem on the unit square.
The PDE is similar to reduced_basis_ex1, except there is a time-derivative
in this case.


<br><br>Basic include file needed for the mesh functionality.
</div>

<div class ="fragment">
<pre>
        #include "libmesh/libmesh.h"
        #include "libmesh/mesh.h"
        #include "libmesh/mesh_generation.h"
        #include "libmesh/exodusII_io.h"
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
        
</pre>
</div>
<div class = "comment">
This example requires SLEPc
</div>

<div class ="fragment">
<pre>
        #if !defined(LIBMESH_HAVE_SLEPC)
          libmesh_example_assert(false, "--enable-slepc");
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
Parse the input file (reduced_basis_ex3.in) using GetPot
</div>

<div class ="fragment">
<pre>
          std::string parameters_filename = "reduced_basis_ex3.in";
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
Finally, we need to give the RBConstruction object a pointer to
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
            unsigned int online_N = infile("online_N",1);
            Real online_x_vel = infile("online_x_vel", 0.);
            Real online_y_vel = infile("online_y_vel", 0.);
            RBParameters online_mu;
            online_mu.set_value("x_vel", online_x_vel);
            online_mu.set_value("y_vel", online_y_vel);
            rb_eval.set_parameters(online_mu);
            rb_eval.print_parameters();
        
</pre>
</div>
<div class = "comment">
Now do the Online solve using the precomputed reduced basis
</div>

<div class ="fragment">
<pre>
            Real error_bound_final_time = rb_eval.rb_solve(online_N);
            
            libMesh::out &lt;&lt; "Error bound (absolute) at the final time is "
                         &lt;&lt; error_bound_final_time &lt;&lt; std::endl &lt;&lt; std::endl;
        
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
Plot the solution at the final time level
</div>

<div class ="fragment">
<pre>
              rb_con.pull_temporal_discretization_data( rb_eval );
              rb_con.set_time_step(rb_con.get_n_time_steps());
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
        
        #endif // LIBMESH_HAVE_SLEPC
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
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/transient_rb_theta_expansion.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/transient_rb_assembly_expansion.h&quot;</FONT></B>
  
  using libMesh::AutoPtr;
  using libMesh::DirichletBoundary;
  using libMesh::ElemAssembly;
  using libMesh::FEInterface;
  using libMesh::FEMContext;
  using libMesh::Number;
  using libMesh::Point;
  using libMesh::RBParameters;
  using libMesh::RBTheta;
  using libMesh::Real;
  using libMesh::RealGradient;
  using libMesh::TransientRBThetaExpansion;
  using libMesh::TransientRBAssemblyExpansion;
  
  <B><FONT COLOR="#228B22">struct</FONT></B> ThetaA0 : RBTheta { <B><FONT COLOR="#228B22">virtual</FONT></B> Number evaluate(<B><FONT COLOR="#228B22">const</FONT></B> RBParameters&amp; )   { <B><FONT COLOR="#A020F0">return</FONT></B> 0.05;  } };
  <B><FONT COLOR="#228B22">struct</FONT></B> ThetaA1 : RBTheta { <B><FONT COLOR="#228B22">virtual</FONT></B> Number evaluate(<B><FONT COLOR="#228B22">const</FONT></B> RBParameters&amp; mu) { <B><FONT COLOR="#A020F0">return</FONT></B> mu.get_value(<B><FONT COLOR="#BC8F8F">&quot;x_vel&quot;</FONT></B>); } };
  <B><FONT COLOR="#228B22">struct</FONT></B> ThetaA2 : RBTheta { <B><FONT COLOR="#228B22">virtual</FONT></B> Number evaluate(<B><FONT COLOR="#228B22">const</FONT></B> RBParameters&amp; mu) { <B><FONT COLOR="#A020F0">return</FONT></B> mu.get_value(<B><FONT COLOR="#BC8F8F">&quot;y_vel&quot;</FONT></B>); } };
  
  <B><FONT COLOR="#228B22">struct</FONT></B> M0 : ElemAssembly
  {
    <B><FONT COLOR="#228B22">virtual</FONT></B> <B><FONT COLOR="#228B22">void</FONT></B> interior_assembly(FEMContext &amp;c)
    {
      <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> u_var = 0;
  
      <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Real&gt; &amp;JxW =
        c.element_fe_var[u_var]-&gt;get_JxW();
  
      <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;Real&gt; &gt;&amp; phi =
        c.element_fe_var[u_var]-&gt;get_phi();
  
      <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_u_dofs = c.dof_indices_var[u_var].size();
  
      <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_qpoints = (c.get_element_qrule())-&gt;n_points();
  
      <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> qp=0; qp != n_qpoints; qp++)
        <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i != n_u_dofs; i++)
          <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j=0; j != n_u_dofs; j++)
            c.get_elem_jacobian()(i,j) += JxW[qp] *phi[j][qp]*phi[i][qp];
    }
  };
  
  <B><FONT COLOR="#228B22">struct</FONT></B> A0 : ElemAssembly
  {
    <B><FONT COLOR="#228B22">virtual</FONT></B> <B><FONT COLOR="#228B22">void</FONT></B> interior_assembly(FEMContext &amp;c)
    {
      <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> u_var = 0;
  
      <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Real&gt; &amp;JxW =
        c.element_fe_var[u_var]-&gt;get_JxW();
  
      <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;RealGradient&gt; &gt;&amp; dphi =
        c.element_fe_var[u_var]-&gt;get_dphi();
  
      <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_u_dofs = c.dof_indices_var[u_var].size();
  
      <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_qpoints = (c.get_element_qrule())-&gt;n_points();
  
      <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> qp=0; qp != n_qpoints; qp++)
        <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i != n_u_dofs; i++)
          <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j=0; j != n_u_dofs; j++)
            c.get_elem_jacobian()(i,j) += JxW[qp] * dphi[j][qp]*dphi[i][qp];
    }
  };
  
  
  <B><FONT COLOR="#228B22">struct</FONT></B> A1 : ElemAssembly
  {
    <B><FONT COLOR="#228B22">virtual</FONT></B> <B><FONT COLOR="#228B22">void</FONT></B> interior_assembly(FEMContext &amp;c)
    {
      <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> u_var = 0;
  
      <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Real&gt; &amp;JxW =
        c.element_fe_var[u_var]-&gt;get_JxW();
  
      <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;Real&gt; &gt;&amp; phi =
        c.element_fe_var[u_var]-&gt;get_phi();
  
      <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;RealGradient&gt; &gt;&amp; dphi =
        c.element_fe_var[u_var]-&gt;get_dphi();
  
      <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_u_dofs = c.dof_indices_var[u_var].size();
  
      <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_qpoints = (c.get_element_qrule())-&gt;n_points();
  
      <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> qp=0; qp != n_qpoints; qp++)
        <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i != n_u_dofs; i++)
          <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j=0; j != n_u_dofs; j++)
            c.get_elem_jacobian()(i,j) += JxW[qp] *dphi[j][qp](0)*phi[i][qp];
    }
  };
  
  <B><FONT COLOR="#228B22">struct</FONT></B> A2 : ElemAssembly
  {
    <B><FONT COLOR="#228B22">virtual</FONT></B> <B><FONT COLOR="#228B22">void</FONT></B> interior_assembly(FEMContext &amp;c)
    {
      <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> u_var = 0;
  
      <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Real&gt; &amp;JxW =
        c.element_fe_var[u_var]-&gt;get_JxW();
  
      <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;Real&gt; &gt;&amp; phi =
        c.element_fe_var[u_var]-&gt;get_phi();
  
      <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;RealGradient&gt; &gt;&amp; dphi =
        c.element_fe_var[u_var]-&gt;get_dphi();
  
      <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_u_dofs = c.dof_indices_var[u_var].size();
  
      <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_qpoints = (c.get_element_qrule())-&gt;n_points();
  
      <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> qp=0; qp != n_qpoints; qp++)
        <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i != n_u_dofs; i++)
          <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j=0; j != n_u_dofs; j++)
            c.get_elem_jacobian()(i,j) += JxW[qp] *dphi[j][qp](1)*phi[i][qp];
    }
  };
  
  <B><FONT COLOR="#228B22">struct</FONT></B> F0 : ElemAssembly
  {
    <B><FONT COLOR="#228B22">virtual</FONT></B> <B><FONT COLOR="#228B22">void</FONT></B> interior_assembly(FEMContext &amp;c)
    {
      <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> u_var = 0;
  
      <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Real&gt; &amp;JxW =
        c.element_fe_var[u_var]-&gt;get_JxW();
  
      <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;Real&gt; &gt;&amp; phi =
        c.element_fe_var[u_var]-&gt;get_phi();
  
      <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_u_dofs = c.dof_indices_var[u_var].size();
  
      <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_qpoints = (c.get_element_qrule())-&gt;n_points();
  
      <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> qp=0; qp != n_qpoints; qp++)
        <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i != n_u_dofs; i++)
          c.get_elem_residual()(i) += JxW[qp] * ( 1.*phi[i][qp] );
    }
  };
  
  <B><FONT COLOR="#228B22">struct</FONT></B> OutputAssembly : ElemAssembly
  {
    OutputAssembly(Real min_x_in, Real max_x_in,
                   Real min_y_in, Real max_y_in)
                  :
                  min_x(min_x_in),
                  max_x(max_x_in),
                  min_y(min_y_in),
                  max_y(max_y_in)
    {}
  
    <B><FONT COLOR="#228B22">virtual</FONT></B> <B><FONT COLOR="#228B22">void</FONT></B> interior_assembly(FEMContext &amp;c)
    {
      <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> u_var = 0;
  
      <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Real&gt; &amp;JxW =
        c.element_fe_var[u_var]-&gt;get_JxW();
  
      <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;Real&gt; &gt;&amp; phi =
        c.element_fe_var[u_var]-&gt;get_phi();
  
      <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_u_dofs = c.dof_indices_var[u_var].size();
  
      <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_qpoints = (c.get_element_qrule())-&gt;n_points();
      
      Real output_area = (max_x-min_x) * (max_y-min_y);
  
      Point centroid = c.elem-&gt;centroid();
      <B><FONT COLOR="#A020F0">if</FONT></B>( (min_x &lt;= centroid(0)) &amp;&amp; (centroid(0) &lt;= max_x) &amp;&amp;
          (min_y &lt;= centroid(1)) &amp;&amp; (centroid(1) &lt;= max_y) )
        <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> qp=0; qp != n_qpoints; qp++)
          <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i != n_u_dofs; i++)
            c.get_elem_residual()(i) += JxW[qp] * ( 1.*phi[i][qp] ) / output_area;
    }
    
    Real min_x, max_x, min_y, max_y;
  };
  
  <B><FONT COLOR="#228B22">struct</FONT></B> CDRBThetaExpansion : TransientRBThetaExpansion
  {
  
    <I><FONT COLOR="#B22222">/**
     * Constructor.
     */</FONT></I>
    CDRBThetaExpansion()
    {
      attach_M_theta(&amp;rb_theta);    <I><FONT COLOR="#B22222">// Attach the time-derivative theta
</FONT></I>  
      attach_A_theta(&amp;theta_a_0);   <I><FONT COLOR="#B22222">// Attach the lhs theta
</FONT></I>      attach_A_theta(&amp;theta_a_1);
      attach_A_theta(&amp;theta_a_2);
  
      attach_F_theta(&amp;rb_theta);    <I><FONT COLOR="#B22222">// Attach the rhs theta
</FONT></I>  
      attach_output_theta(&amp;rb_theta); <I><FONT COLOR="#B22222">// Attach output 0 theta
</FONT></I>      attach_output_theta(&amp;rb_theta); <I><FONT COLOR="#B22222">// Attach output 1 theta
</FONT></I>      attach_output_theta(&amp;rb_theta); <I><FONT COLOR="#B22222">// Attach output 2 theta
</FONT></I>      attach_output_theta(&amp;rb_theta); <I><FONT COLOR="#B22222">// Attach output 3 theta
</FONT></I>    }
  
    ThetaA0 theta_a_0;
    ThetaA1 theta_a_1;
    ThetaA2 theta_a_2;
    RBTheta rb_theta; <I><FONT COLOR="#B22222">// Default RBTheta object, just returns 1.
</FONT></I>  };
  
  <B><FONT COLOR="#228B22">struct</FONT></B> CDRBAssemblyExpansion : TransientRBAssemblyExpansion
  {
  
    <I><FONT COLOR="#B22222">/**
     * Constructor.
     */</FONT></I>
    CDRBAssemblyExpansion()
      :
      L0(0.72,0.88,0.72,0.88), <I><FONT COLOR="#B22222">// We make sure these output regions conform to the mesh
</FONT></I>      L1(0.12,0.28,0.72,0.88),
      L2(0.12,0.28,0.12,0.28),
      L3(0.72,0.88,0.12,0.28)
    {
      attach_M_assembly(&amp;M0_assembly); <I><FONT COLOR="#B22222">// Attach the time-derivative assembly
</FONT></I>      attach_A_assembly(&amp;A0_assembly); <I><FONT COLOR="#B22222">// Attach the lhs assembly
</FONT></I>      attach_A_assembly(&amp;A1_assembly);
      attach_A_assembly(&amp;A2_assembly);
      
      attach_F_assembly(&amp;F0_assembly); <I><FONT COLOR="#B22222">// Attach the rhs assembly
</FONT></I>      
      attach_output_assembly(&amp;L0);       <I><FONT COLOR="#B22222">// Attach output 0 assembly
</FONT></I>      attach_output_assembly(&amp;L1);       <I><FONT COLOR="#B22222">// Attach output 1 assembly
</FONT></I>      attach_output_assembly(&amp;L2);       <I><FONT COLOR="#B22222">// Attach output 2 assembly
</FONT></I>      attach_output_assembly(&amp;L3);       <I><FONT COLOR="#B22222">// Attach output 3 assembly
</FONT></I>    }
  
    M0 M0_assembly;
    A0 A0_assembly;
    A1 A1_assembly;
    A2 A2_assembly;
    F0 F0_assembly;
    OutputAssembly L0;
    OutputAssembly L1;
    OutputAssembly L2;
    OutputAssembly L3;
  };
  
  #endif
</pre> 
<a name="nocomments"></a> 
<br><br><br> <h1> The source file rb_classes.h without comments: </h1> 
<pre> 
    
    
  
  #ifndef __rb_classes_h__
  #define __rb_classes_h__
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/transient_rb_construction.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/fe_base.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;assembly.h&quot;</FONT></B>
  
  using libMesh::EquationSystems;
  using libMesh::FEMContext;
  using libMesh::RBConstruction;
  using libMesh::RBEvaluation;
  using libMesh::Real;
  using libMesh::TransientRBEvaluation;
  using libMesh::TransientRBConstruction;
  
  
  <B><FONT COLOR="#228B22">class</FONT></B> SimpleRBEvaluation : <B><FONT COLOR="#228B22">public</FONT></B> TransientRBEvaluation
  {
  <B><FONT COLOR="#228B22">public</FONT></B>:
  
    <I><FONT COLOR="#B22222">/**
     * Constructor. Just set the theta expansion.
     */</FONT></I>
    SimpleRBEvaluation()
    {
      set_rb_theta_expansion(cd_rb_theta_expansion);
    }
  
    <I><FONT COLOR="#B22222">/**
     * The coercivity constant is bounded below by 0.05.
     */</FONT></I>
    <B><FONT COLOR="#228B22">virtual</FONT></B> Real get_stability_lower_bound() { <B><FONT COLOR="#A020F0">return</FONT></B> 0.05; }
  
    <I><FONT COLOR="#B22222">/**
     * The object that stores the &quot;theta&quot; expansion of the parameter dependent PDE,
     * i.e. the set of parameter-dependent functions in the affine expansion of the PDE.
     */</FONT></I>
    CDRBThetaExpansion cd_rb_theta_expansion;
  
  };
  
  <B><FONT COLOR="#228B22">class</FONT></B> SimpleRBConstruction : <B><FONT COLOR="#228B22">public</FONT></B> TransientRBConstruction
  {
  <B><FONT COLOR="#228B22">public</FONT></B>:
  
    SimpleRBConstruction (EquationSystems&amp; es,
                          <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp; name,
                          <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> number)
    : Parent(es, name, number),
      dirichlet_bc(AutoPtr&lt;DirichletBoundary&gt;(NULL))
    {}
  
    <I><FONT COLOR="#B22222">/**
     * Destructor.
     */</FONT></I>
    <B><FONT COLOR="#228B22">virtual</FONT></B> ~SimpleRBConstruction () {}
  
    <I><FONT COLOR="#B22222">/**
     * The type of system.
     */</FONT></I>
    <B><FONT COLOR="#228B22">typedef</FONT></B> SimpleRBConstruction sys_type;
  
    <I><FONT COLOR="#B22222">/**
     * The type of the parent.
     */</FONT></I>
    <B><FONT COLOR="#228B22">typedef</FONT></B> TransientRBConstruction Parent;
  
    <I><FONT COLOR="#B22222">/**
     * Initialize data structures.
     */</FONT></I>
    <B><FONT COLOR="#228B22">virtual</FONT></B> <B><FONT COLOR="#228B22">void</FONT></B> init_data()
    {
      u_var = <B><FONT COLOR="#A020F0">this</FONT></B>-&gt;add_variable (<B><FONT COLOR="#BC8F8F">&quot;u&quot;</FONT></B>, FIRST);
  
      dirichlet_bc = build_zero_dirichlet_boundary_object();
  
      dirichlet_bc-&gt;b.insert(0);
      dirichlet_bc-&gt;b.insert(1);
      dirichlet_bc-&gt;b.insert(2);
      dirichlet_bc-&gt;b.insert(3);
      dirichlet_bc-&gt;variables.push_back(u_var);
      
      get_dof_map().add_dirichlet_boundary(*dirichlet_bc);
  
      <B><FONT COLOR="#5F9EA0">Parent</FONT></B>::init_data();
  
      set_rb_assembly_expansion(cd_rb_assembly_expansion);
  
      set_inner_product_assembly(cd_rb_assembly_expansion.A0_assembly);
  
      set_L2_assembly(cd_rb_assembly_expansion.M0_assembly);
    }
  
    <I><FONT COLOR="#B22222">/**
     * Pre-request all relevant element data.
     */</FONT></I>
    <B><FONT COLOR="#228B22">virtual</FONT></B> <B><FONT COLOR="#228B22">void</FONT></B> init_context(FEMContext &amp;c)
    {
      c.element_fe_var[u_var]-&gt;get_JxW();
      c.element_fe_var[u_var]-&gt;get_phi();
      c.element_fe_var[u_var]-&gt;get_dphi();
    }
  
    <I><FONT COLOR="#B22222">/**
     * Variable number for u.
     */</FONT></I>
    <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> u_var;
    
    <I><FONT COLOR="#B22222">/**
     * The object that stores the &quot;assembly&quot; expansion of the parameter dependent PDE,
     * i.e. the objects that define how to assemble the set of parameter-independent
     * operators in the affine expansion of the PDE.
     */</FONT></I>
    CDRBAssemblyExpansion cd_rb_assembly_expansion;
  
    <I><FONT COLOR="#B22222">/**
     * The object that defines which degrees of freedom are on a Dirichlet boundary.
     */</FONT></I>
    AutoPtr&lt;DirichletBoundary&gt; dirichlet_bc;
  
  };
  
  #endif
</pre> 
<a name="nocomments"></a> 
<br><br><br> <h1> The source file reduced_basis_ex3.C without comments: </h1> 
<pre> 
    
  <I><FONT COLOR="#B22222">/* rbOOmit is distributed in the hope that it will be useful, */</FONT></I>
  <I><FONT COLOR="#B22222">/* but WITHOUT ANY WARRANTY; without even the implied warranty of */</FONT></I>
  <I><FONT COLOR="#B22222">/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU */</FONT></I>
  <I><FONT COLOR="#B22222">/* Lesser General Public License for more details. */</FONT></I>
    
  <I><FONT COLOR="#B22222">/* You should have received a copy of the GNU Lesser General Public */</FONT></I>
  <I><FONT COLOR="#B22222">/* License along with this library; if not, write to the Free Software */</FONT></I>
  <I><FONT COLOR="#B22222">/* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA */</FONT></I>
  
  
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/libmesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/mesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/mesh_generation.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/exodusII_io.h&quot;</FONT></B>
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
  
  #<B><FONT COLOR="#A020F0">if</FONT></B> !defined(LIBMESH_HAVE_SLEPC)
    libmesh_example_assert(false, <B><FONT COLOR="#BC8F8F">&quot;--enable-slepc&quot;</FONT></B>);
  #<B><FONT COLOR="#A020F0">else</FONT></B>
  
  #<B><FONT COLOR="#A020F0">if</FONT></B> !defined(LIBMESH_HAVE_XDR)
    libmesh_example_assert(false, <B><FONT COLOR="#BC8F8F">&quot;--enable-xdr&quot;</FONT></B>);
  #elif defined(LIBMESH_DEFAULT_SINGLE_PRECISION)
    libmesh_example_assert(false, <B><FONT COLOR="#BC8F8F">&quot;--disable-singleprecision&quot;</FONT></B>);
  #endif
    libmesh_example_assert(libMesh::default_solver_package() == PETSC_SOLVERS, <B><FONT COLOR="#BC8F8F">&quot;--enable-petsc&quot;</FONT></B>);
  
    libmesh_example_assert(2 &lt;= LIBMESH_DIM, <B><FONT COLOR="#BC8F8F">&quot;2D support&quot;</FONT></B>);
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::string parameters_filename = <B><FONT COLOR="#BC8F8F">&quot;reduced_basis_ex3.in&quot;</FONT></B>;
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
      
      <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> online_N = infile(<B><FONT COLOR="#BC8F8F">&quot;online_N&quot;</FONT></B>,1);
      Real online_x_vel = infile(<B><FONT COLOR="#BC8F8F">&quot;online_x_vel&quot;</FONT></B>, 0.);
      Real online_y_vel = infile(<B><FONT COLOR="#BC8F8F">&quot;online_y_vel&quot;</FONT></B>, 0.);
      RBParameters online_mu;
      online_mu.set_value(<B><FONT COLOR="#BC8F8F">&quot;x_vel&quot;</FONT></B>, online_x_vel);
      online_mu.set_value(<B><FONT COLOR="#BC8F8F">&quot;y_vel&quot;</FONT></B>, online_y_vel);
      rb_eval.set_parameters(online_mu);
      rb_eval.print_parameters();
  
      Real error_bound_final_time = rb_eval.rb_solve(online_N);
      
      <B><FONT COLOR="#5F9EA0">libMesh</FONT></B>::out &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Error bound (absolute) at the final time is &quot;</FONT></B>
                   &lt;&lt; error_bound_final_time &lt;&lt; std::endl &lt;&lt; std::endl;
  
      <B><FONT COLOR="#A020F0">if</FONT></B>(store_basis_functions)
      {
        rb_eval.read_in_basis_functions(rb_con);
        
        rb_con.pull_temporal_discretization_data( rb_eval );
        rb_con.set_time_step(rb_con.get_n_time_steps());
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
  
  #endif <I><FONT COLOR="#B22222">// LIBMESH_HAVE_SLEPC
</FONT></I>  }
  
</pre> 
<a name="output"></a> 
<br><br><br> <h1> The console output of the program: </h1> 
<pre>
***************************************************************
* Running Example reduced_basis_ex3:
*  mpirun -np 12 example-devel -online_mode 0 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
***************************************************************
 
*** Warning, This code is untested, experimental, or likely to see future API changes: src/reduced_basis/rb_parametrized.C, line 41, compiled Jan 31 2013 at 21:51:32 ***
 EquationSystems
  n_systems()=1
   System #0, "RBConvectionDiffusion"
    Type "TransientRBConstruction"
    Variables="u" 
    Finite Element Types="LAGRANGE", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" 
    Approximation Orders="FIRST", "THIRD" 
    n_dofs()=676
    n_local_dofs()=68
    n_constrained_dofs()=100
    n_local_constrained_dofs()=16
    n_vectors()=3
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 7.99852
      Average Off-Processor Bandwidth <= 1.20118
      Maximum  On-Processor Bandwidth <= 11
      Maximum Off-Processor Bandwidth <= 9
    DofMap Constraints
      Number of DoF Constraints = 100
      Average DoF Constraint Length= 0
      Number of Node Constraints = 0

 Mesh Information:
  mesh_dimension()=2
  spatial_dimension()=3
  n_nodes()=676
    n_local_nodes()=68
  n_elem()=625
    n_local_elem()=51
    n_active_elem()=625
  n_subdomains()=1
  n_partitions()=12
  n_processors()=12
  n_threads()=1
  processor_id()=0

Initializing training parameters with deterministic training set...
Parameter x_vel: log scaling = 0
Parameter y_vel: log scaling = 0


RBConstruction parameters:
system name: RBConvectionDiffusion
constrained_problem: 0
Nmax: 20
Aq operators attached: 3
Fq functions attached: 1
n_outputs: 4
output 0, Q_l = 1
output 1, Q_l = 1
output 2, Q_l = 1
output 3, Q_l = 1
Number of parameters: 2
Parameter x_vel: Min = -2, Max = 2, value = 1
Parameter y_vel: Min = -2, Max = 2, value = 1
n_training_samples: 100
single-matrix mode? 0
reuse preconditioner? 1
use a relative error bound in greedy? 1
write out data during basis training? 0
quiet mode? 1


TransientRBConstruction parameters:
Q_m: 1
Number of time-steps: 100
dt: 0.01
euler_theta (time discretization parameter): 1
delta_N (number of basis functions to add each POD-Greedy step): 1
Using zero initial condition

Compute output dual inner products
output_dual_innerprods[0][0] = 33.6538
output_dual_innerprods[1][0] = 33.6538
output_dual_innerprods[2][0] = 33.6538
output_dual_innerprods[3][0] = 33.6538

---- Performing Greedy basis enrichment ----

---- Basis dimension: 0 ----
Performing RB solves on training set
Maximum (relative) error bound is inf

Performing truth solve at parameter:
x_vel: -2.000000e+00
y_vel: -2.000000e+00

Enriching the RB space

POD Eigenvalues:
eigenvalue 0 = 2.19193
eigenvalue 1 = 0.0116709
eigenvalue 2 = 0.00106473
...
last eigenvalue = -1.19071e-16

Updating RB matrices
Updating RB residual terms
Updating RB initial conditions

---- Basis dimension: 1 ----
Performing RB solves on training set
Maximum (relative) error bound is 6.81811

Performing truth solve at parameter:
x_vel: 2.000000e+00
y_vel: 2.000000e+00

Enriching the RB space

POD Eigenvalues:
eigenvalue 0 = 2.14617
eigenvalue 1 = 0.0106903
eigenvalue 2 = 0.000908038
...
last eigenvalue = -1.93388e-16

Updating RB matrices
Updating RB residual terms
Updating RB initial conditions

---- Basis dimension: 2 ----
Performing RB solves on training set
Maximum (relative) error bound is 3.15659

Performing truth solve at parameter:
x_vel: 2.000000e+00
y_vel: -2.000000e+00

Enriching the RB space

POD Eigenvalues:
eigenvalue 0 = 1.6024
eigenvalue 1 = 0.00595727
eigenvalue 2 = 0.000820407
...
last eigenvalue = -9.72915e-17

Updating RB matrices
Updating RB residual terms
Updating RB initial conditions

---- Basis dimension: 3 ----
Performing RB solves on training set
Maximum (relative) error bound is 3.95604

Performing truth solve at parameter:
x_vel: -2.000000e+00
y_vel: 2.000000e+00

Enriching the RB space

POD Eigenvalues:
eigenvalue 0 = 1.55558
eigenvalue 1 = 0.0058283
eigenvalue 2 = 0.00079201
...
last eigenvalue = -9.04825e-17

Updating RB matrices
Updating RB residual terms
Updating RB initial conditions

---- Basis dimension: 4 ----
Performing RB solves on training set
Maximum (relative) error bound is 1.3165

Performing truth solve at parameter:
x_vel: -2.222222e-01
y_vel: 2.222222e-01

Enriching the RB space

POD Eigenvalues:
eigenvalue 0 = 1.02449
eigenvalue 1 = 0.0148056
eigenvalue 2 = 0.000593638
...
last eigenvalue = -1.61741e-16

Updating RB matrices
Updating RB residual terms
Updating RB initial conditions

---- Basis dimension: 5 ----
Performing RB solves on training set
Maximum (relative) error bound is 1.05097

Performing truth solve at parameter:
x_vel: -2.000000e+00
y_vel: -2.222222e-01

Enriching the RB space

POD Eigenvalues:
eigenvalue 0 = 0.365693
eigenvalue 1 = 0.00367148
eigenvalue 2 = 0.000675947
...
last eigenvalue = -6.01347e-18

Updating RB matrices
Updating RB residual terms
Updating RB initial conditions

---- Basis dimension: 6 ----
Performing RB solves on training set
Maximum (relative) error bound is 1.0171

Performing truth solve at parameter:
x_vel: -2.222222e-01
y_vel: -2.000000e+00

Enriching the RB space

POD Eigenvalues:
eigenvalue 0 = 0.285715
eigenvalue 1 = 0.00305589
eigenvalue 2 = 0.000804404
...
last eigenvalue = -6.23967e-17

Updating RB matrices
Updating RB residual terms
Updating RB initial conditions

---- Basis dimension: 7 ----
Performing RB solves on training set
Maximum (relative) error bound is 1.14302

Performing truth solve at parameter:
x_vel: -2.222222e-01
y_vel: 2.000000e+00

Enriching the RB space

POD Eigenvalues:
eigenvalue 0 = 0.256366
eigenvalue 1 = 0.00345875
eigenvalue 2 = 0.000587394
...
last eigenvalue = -2.28299e-18

Updating RB matrices
Updating RB residual terms
Updating RB initial conditions

---- Basis dimension: 8 ----
Performing RB solves on training set
Maximum (relative) error bound is 1.21494

Performing truth solve at parameter:
x_vel: 2.000000e+00
y_vel: -2.222222e-01

Enriching the RB space

POD Eigenvalues:
eigenvalue 0 = 0.260659
eigenvalue 1 = 0.00277678
eigenvalue 2 = 0.000773864
...
last eigenvalue = -1.6619e-17

Updating RB matrices
Updating RB residual terms
Updating RB initial conditions

---- Basis dimension: 9 ----
Performing RB solves on training set
Maximum (relative) error bound is 0.43786

Performing truth solve at parameter:
x_vel: 1.111111e+00
y_vel: 1.111111e+00

Enriching the RB space

POD Eigenvalues:
eigenvalue 0 = 0.0933151
eigenvalue 1 = 0.00402501
eigenvalue 2 = 0.00138016
...
last eigenvalue = -1.01915e-17

Updating RB matrices
Updating RB residual terms
Updating RB initial conditions

---- Basis dimension: 10 ----
Performing RB solves on training set
Maximum (relative) error bound is 0.433787

Performing truth solve at parameter:
x_vel: 6.666667e-01
y_vel: -1.111111e+00

Enriching the RB space

POD Eigenvalues:
eigenvalue 0 = 0.134913
eigenvalue 1 = 0.00436435
eigenvalue 2 = 0.00137045
...
last eigenvalue = -1.72368e-17

Updating RB matrices
Updating RB residual terms
Updating RB initial conditions

---- Basis dimension: 11 ----
Performing RB solves on training set
Maximum (relative) error bound is 0.319326

Performing truth solve at parameter:
x_vel: -1.111111e+00
y_vel: -1.111111e+00

Enriching the RB space

POD Eigenvalues:
eigenvalue 0 = 0.0629291
eigenvalue 1 = 0.00315345
eigenvalue 2 = 0.00120446
...
last eigenvalue = -6.81126e-18

Updating RB matrices
Updating RB residual terms
Updating RB initial conditions

---- Basis dimension: 12 ----
Performing RB solves on training set
Maximum (relative) error bound is 0.287715

Performing truth solve at parameter:
x_vel: -1.555556e+00
y_vel: 1.111111e+00

Enriching the RB space

POD Eigenvalues:
eigenvalue 0 = 0.0278364
eigenvalue 1 = 0.00277066
eigenvalue 2 = 0.000760673
...
last eigenvalue = -2.86335e-19

Updating RB matrices
Updating RB residual terms
Updating RB initial conditions

---- Basis dimension: 13 ----
Performing RB solves on training set
Maximum (relative) error bound is 0.228017

Performing truth solve at parameter:
x_vel: 1.111111e+00
y_vel: -2.000000e+00

Enriching the RB space

POD Eigenvalues:
eigenvalue 0 = 0.00863589
eigenvalue 1 = 0.00176423
eigenvalue 2 = 0.000525486
...
last eigenvalue = -1.87818e-18

Updating RB matrices
Updating RB residual terms
Updating RB initial conditions

---- Basis dimension: 14 ----
Performing RB solves on training set
Maximum (relative) error bound is 0.227118

Performing truth solve at parameter:
x_vel: 2.000000e+00
y_vel: 1.111111e+00

Enriching the RB space

POD Eigenvalues:
eigenvalue 0 = 0.0117352
eigenvalue 1 = 0.00206357
eigenvalue 2 = 0.000546886
...
last eigenvalue = -9.40987e-19

Updating RB matrices
Updating RB residual terms
Updating RB initial conditions

---- Basis dimension: 15 ----
Performing RB solves on training set
Maximum (relative) error bound is 0.231987

Performing truth solve at parameter:
x_vel: 6.666667e-01
y_vel: 2.000000e+00

Enriching the RB space

POD Eigenvalues:
eigenvalue 0 = 0.00898933
eigenvalue 1 = 0.0022199
eigenvalue 2 = 0.000417343
...
last eigenvalue = -8.9785e-19

Updating RB matrices
Updating RB residual terms
Updating RB initial conditions

---- Basis dimension: 16 ----
Performing RB solves on training set
Maximum (relative) error bound is 0.199807

Performing truth solve at parameter:
x_vel: -2.000000e+00
y_vel: 6.666667e-01

Enriching the RB space

POD Eigenvalues:
eigenvalue 0 = 0.0065467
eigenvalue 1 = 0.00176326
eigenvalue 2 = 0.000231436
...
last eigenvalue = -4.2787e-19

Updating RB matrices
Updating RB residual terms
Updating RB initial conditions

---- Basis dimension: 17 ----
Performing RB solves on training set
Maximum (relative) error bound is 0.193551

Performing truth solve at parameter:
x_vel: -1.111111e+00
y_vel: 2.000000e+00

Enriching the RB space

POD Eigenvalues:
eigenvalue 0 = 0.00390004
eigenvalue 1 = 0.00190681
eigenvalue 2 = 0.000456244
...
last eigenvalue = -1.71625e-19

Updating RB matrices
Updating RB residual terms
Updating RB initial conditions

---- Basis dimension: 18 ----
Performing RB solves on training set
Maximum (relative) error bound is 0.181413

Performing truth solve at parameter:
x_vel: -1.111111e+00
y_vel: -2.000000e+00

Enriching the RB space

POD Eigenvalues:
eigenvalue 0 = 0.00391101
eigenvalue 1 = 0.00171004
eigenvalue 2 = 0.000346167
...
last eigenvalue = -9.60514e-20

Updating RB matrices
Updating RB residual terms
Updating RB initial conditions

---- Basis dimension: 19 ----
Performing RB solves on training set
Maximum (relative) error bound is 0.174248

Performing truth solve at parameter:
x_vel: 1.111111e+00
y_vel: -2.222222e-01

Enriching the RB space

POD Eigenvalues:
eigenvalue 0 = 0.0287899
eigenvalue 1 = 0.00198823
eigenvalue 2 = 0.000531507
...
last eigenvalue = -8.67362e-19

Updating RB matrices
Updating RB residual terms
Updating RB initial conditions

---- Basis dimension: 20 ----
Performing RB solves on training set
Maximum (relative) error bound is 0.145574

Maximum number of basis functions reached: Nmax = 20
************************************************************************************************************************
***             WIDEN YOUR WINDOW TO 120 CHARACTERS.  Use 'enscript -r -fCourier9' to print this document            ***
************************************************************************************************************************

---------------------------------------------- PETSc Performance Summary: ----------------------------------------------

/workspace/libmesh/examples/reduced_basis/reduced_basis_ex3/.libs/lt-example-devel on a intel-12. named hbar.ices.utexas.edu with 12 processors, by benkirk Thu Jan 31 22:17:34 2013
Using Petsc Release Version 3.3.0, Patch 2, Fri Jul 13 15:42:00 CDT 2012 

                         Max       Max/Min        Avg      Total 
Time (sec):           1.270e+01      1.00000   1.270e+01
Objects:              7.293e+04      1.00000   7.293e+04
Flops:                1.946e+08      1.86315   1.558e+08  1.869e+09
Flops/sec:            1.532e+07      1.86315   1.226e+07  1.472e+08
MPI Messages:         4.033e+05      3.01546   2.682e+05  3.218e+06
MPI Message Lengths:  1.838e+07      2.28703   4.649e+01  1.496e+08
MPI Reductions:       4.262e+05      1.00000

Flop counting convention: 1 flop = 1 real number operation of type (multiply/divide/add/subtract)
                            e.g., VecAXPY() for real vectors of length N --> 2N flops
                            and VecAXPY() for complex vectors of length N --> 8N flops

Summary of Stages:   ----- Time ------  ----- Flops -----  --- Messages ---  -- Message Lengths --  -- Reductions --
                        Avg     %Total     Avg     %Total   counts   %Total     Avg         %Total   counts   %Total 
 0:      Main Stage: 1.2701e+01 100.0%  1.8691e+09 100.0%  3.218e+06 100.0%  4.649e+01      100.0%  4.262e+05 100.0% 

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

VecDot            138405 1.0 8.2195e-01 1.1 1.87e+07 1.8 0.0e+00 0.0e+00 1.4e+05  6 10  0  0 32   6 10  0  0 32   226
VecMDot            24464 1.0 3.4280e-01 1.5 2.42e+07 1.8 0.0e+00 0.0e+00 2.4e+04  2 13  0  0  6   2 13  0  0  6   700
VecNorm            28715 1.0 3.0868e-01 1.2 3.91e+06 1.8 0.0e+00 0.0e+00 2.9e+04  2  2  0  0  7   2  2  0  0  7   126
VecScale           34710 1.0 1.8898e-02 1.3 2.22e+06 1.8 0.0e+00 0.0e+00 0.0e+00  0  1  0  0  0   0  1  0  0  0  1170
VecCopy            10295 1.0 9.5944e-03 1.4 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecSet             73266 1.0 2.9309e-02 1.3 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecAXPY            37797 1.0 3.1516e-02 1.5 5.08e+06 1.8 0.0e+00 0.0e+00 0.0e+00  0  3  0  0  0   0  3  0  0  0  1603
VecMAXPY           26630 1.0 2.4882e-02 1.4 2.77e+07 1.8 0.0e+00 0.0e+00 0.0e+00  0 15  0  0  0   0 15  0  0  0 11062
VecAssemblyBegin   14344 1.0 5.8009e-01 1.2 0.00e+00 0.0 2.4e+02 8.1e+01 4.3e+04  4  0  0  0 10   4  0  0  0 10     0
VecAssemblyEnd     14344 1.0 8.9967e-03 1.3 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecScatterBegin    50749 1.0 1.9577e-01 1.6 0.00e+00 0.0 2.4e+06 5.6e+01 0.0e+00  1  0 76 92  0   1  0 76 92  0     0
VecScatterEnd      50749 1.0 3.3280e-01 2.4 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  2  0  0  0  0   2  0  0  0  0     0
VecNormalize       26630 1.0 3.0794e-01 1.2 5.43e+06 1.8 0.0e+00 0.0e+00 2.7e+04  2  3  0  0  6   2  3  0  0  6   175
MatMult            26630 1.0 3.1586e-01 1.7 2.81e+07 1.6 1.3e+06 5.4e+01 0.0e+00  2 15 40 46  0   2 15 40 46  0   917
MatMultAdd         19994 1.0 2.4625e-01 1.6 2.25e+07 1.6 9.6e+05 5.4e+01 0.0e+00  2 12 30 35  0   2 12 30 35  0   938
MatSolve           28715 1.0 9.1805e-02 2.3 6.19e+07 2.2 0.0e+00 0.0e+00 0.0e+00  0 29  0  0  0   0 29  0  0  0  5815
MatLUFactorNum        42 1.0 1.4949e-03 1.7 3.00e+05 2.5 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0  1641
MatILUFactorSym       42 1.0 3.3593e-03 1.7 0.00e+00 0.0 0.0e+00 0.0e+00 1.3e+02  0  0  0  0  0   0  0  0  0  0     0
MatAssemblyBegin   36214 1.0 1.3596e+00 1.4 0.00e+00 0.0 4.3e+02 2.5e+02 7.2e+04  9  0  0  0 17   9  0  0  0 17     0
MatAssemblyEnd     36214 1.0 1.5418e+00 1.1 0.00e+00 0.0 7.7e+05 1.6e+01 6.4e+04 12  0 24  8 15  12  0 24  8 15     0
MatGetRow        1093712 1.8 2.0853e-01 1.6 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  1  0  0  0  0   1  0  0  0  0     0
MatGetRowIJ           42 1.0 2.3365e-05 2.7 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatGetOrdering        42 1.0 1.3680e-03 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 8.4e+01  0  0  0  0  0   0  0  0  0  0     0
MatZeroEntries      2062 1.0 2.1219e-03 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatAXPY             8042 1.0 5.8290e+00 1.0 0.00e+00 0.0 7.7e+05 1.6e+01 1.3e+05 46  0 24  8 30  46  0 24  8 30     0
KSPGMRESOrthog     24464 1.0 3.8382e-01 1.4 4.85e+07 1.8 0.0e+00 0.0e+00 2.4e+04  3 26  0  0  6   3 26  0  0  6  1256
KSPSetUp            2127 1.0 8.9834e-03 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPSolve            2085 1.0 1.3315e+00 1.0 1.49e+08 1.9 1.3e+06 5.4e+01 5.3e+04 10 75 40 46 13  10 75 40 46 13  1054
PCSetUp               84 1.0 1.6003e-02 1.2 3.00e+05 2.5 0.0e+00 0.0e+00 2.1e+02  0  0  0  0  0   0  0  0  0  0   153
PCSetUpOnBlocks     2085 1.0 1.4232e-02 1.2 3.00e+05 2.5 0.0e+00 0.0e+00 2.1e+02  0  0  0  0  0   0  0  0  0  0   172
PCApply            28715 1.0 3.4235e-01 1.2 6.19e+07 2.2 0.0e+00 0.0e+00 0.0e+00  2 29  0  0  0   2 29  0  0  0  1559
------------------------------------------------------------------------------------------------------------------------

Memory usage is given in bytes:

Object Type          Creations   Destructions     Memory  Descendants' Mem.
Reports information only for process 0.

--- Event Stage 0: Main Stage

              Vector 24441          24441     41870472     0
      Vector Scatter  8054           8054      8343944     0
           Index Set 16234          16234     12669328     0
   IS L to G Mapping     5              5         2820     0
              Matrix 24189          24189    133868532     0
       Krylov Solver     2              2        19360     0
      Preconditioner     2              2         1784     0
              Viewer     1              0            0     0
========================================================================================================================
Average time to get PetscTime(): 9.53674e-08
Average time for MPI_Barrier(): 3.76701e-06
Average time for zero size MPI_Send(): 1.39078e-05
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
| Time:           Thu Jan 31 22:17:34 2013                                                                             |
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
 -------------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=12.7284, Active time=12.6475                                                      |
 -------------------------------------------------------------------------------------------------------------------
| Event                                 nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                                 w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|-------------------------------------------------------------------------------------------------------------------|
|                                                                                                                   |
|                                                                                                                   |
| DofMap                                                                                                            |
|   add_neighbors_to_send_list()        1         0.0031      0.003137    0.0051      0.005059    0.02     0.04     |
|   build_constraint_matrix()           561       0.0055      0.000010    0.0055      0.000010    0.04     0.04     |
|   build_sparsity()                    1         0.0021      0.002143    0.0067      0.006735    0.02     0.05     |
|   cnstrn_elem_mat_vec()               561       0.0249      0.000044    0.0249      0.000044    0.20     0.20     |
|   create_dof_constraints()            1         0.0121      0.012127    0.0468      0.046811    0.10     0.37     |
|   distribute_dofs()                   1         0.0112      0.011189    0.0365      0.036460    0.09     0.29     |
|   dof_indices()                       1833      0.1332      0.000073    0.1332      0.000073    1.05     1.05     |
|   prepare_send_list()                 1         0.0000      0.000037    0.0000      0.000037    0.00     0.00     |
|   reinit()                            1         0.0235      0.023474    0.0235      0.023474    0.19     0.19     |
|                                                                                                                   |
| FE                                                                                                                |
|   compute_shape_functions()           1452      0.0376      0.000026    0.0376      0.000026    0.30     0.30     |
|   init_shape_functions()              352       0.0121      0.000035    0.0121      0.000035    0.10     0.10     |
|   inverse_map()                       660       0.0074      0.000011    0.0074      0.000011    0.06     0.06     |
|                                                                                                                   |
| FEMap                                                                                                             |
|   compute_affine_map()                1452      0.0148      0.000010    0.0148      0.000010    0.12     0.12     |
|   compute_face_map()                  330       0.0070      0.000021    0.0146      0.000044    0.06     0.12     |
|   init_face_shape_functions()         22        0.0004      0.000018    0.0004      0.000018    0.00     0.00     |
|   init_reference_to_physical_map()    352       0.0072      0.000020    0.0072      0.000020    0.06     0.06     |
|                                                                                                                   |
| Mesh                                                                                                              |
|   find_neighbors()                    1         0.0139      0.013891    0.0146      0.014569    0.11     0.12     |
|   renumber_nodes_and_elem()           2         0.0006      0.000286    0.0006      0.000286    0.00     0.00     |
|                                                                                                                   |
| MeshCommunication                                                                                                 |
|   assign_global_indices()             1         0.0275      0.027503    0.0298      0.029825    0.22     0.24     |
|   compute_hilbert_indices()           2         0.0107      0.005353    0.0107      0.005353    0.08     0.08     |
|   find_global_indices()               2         0.0042      0.002114    0.0202      0.010098    0.03     0.16     |
|   parallel_sort()                     2         0.0031      0.001534    0.0039      0.001943    0.02     0.03     |
|                                                                                                                   |
| MeshTools::Generation                                                                                             |
|   build_cube()                        1         0.0031      0.003141    0.0031      0.003141    0.02     0.02     |
|                                                                                                                   |
| MetisPartitioner                                                                                                  |
|   partition()                         1         0.0494      0.049351    0.0588      0.058775    0.39     0.46     |
|                                                                                                                   |
| Parallel                                                                                                          |
|   allgather()                         14        0.0012      0.000084    0.0013      0.000090    0.01     0.01     |
|   barrier()                           1         0.0000      0.000020    0.0000      0.000020    0.00     0.00     |
|   broadcast()                         44        0.0003      0.000007    0.0003      0.000007    0.00     0.00     |
|   max(bool)                           7         0.0000      0.000006    0.0000      0.000006    0.00     0.00     |
|   max(scalar)                         127       0.0011      0.000009    0.0011      0.000009    0.01     0.01     |
|   max(vector)                         30        0.0005      0.000016    0.0012      0.000041    0.00     0.01     |
|   maxloc(scalar)                      21        0.1808      0.008608    0.1808      0.008608    1.43     1.43     |
|   min(bool)                           148       0.0012      0.000008    0.0012      0.000008    0.01     0.01     |
|   min(scalar)                         121       0.0104      0.000086    0.0104      0.000086    0.08     0.08     |
|   min(vector)                         30        0.0006      0.000019    0.0014      0.000047    0.00     0.01     |
|   probe()                             272       0.0014      0.000005    0.0014      0.000005    0.01     0.01     |
|   receive()                           224       0.0016      0.000007    0.0028      0.000013    0.01     0.02     |
|   send()                              180       0.0006      0.000003    0.0006      0.000003    0.00     0.00     |
|   send_receive()                      184       0.0015      0.000008    0.0048      0.000026    0.01     0.04     |
|   sum()                               44        0.0028      0.000064    0.0036      0.000083    0.02     0.03     |
|                                                                                                                   |
| Parallel::Request                                                                                                 |
|   wait()                              180       0.0004      0.000002    0.0004      0.000002    0.00     0.00     |
|                                                                                                                   |
| Partitioner                                                                                                       |
|   set_node_processor_ids()            1         0.0016      0.001559    0.0023      0.002344    0.01     0.02     |
|   set_parent_processor_ids()          1         0.0012      0.001190    0.0012      0.001190    0.01     0.01     |
|                                                                                                                   |
| PetscLinearSolver                                                                                                 |
|   solve()                             2085      1.6694      0.000801    1.6694      0.000801    13.20    13.20    |
|                                                                                                                   |
| RBConstruction                                                                                                    |
|   add_scaled_matrix_and_vector()      11        0.1979      0.017988    0.4105      0.037317    1.56     3.25     |
|   clear()                             3         0.0013      0.000426    0.0013      0.000426    0.01     0.01     |
|   compute_Fq_representor_innerprods() 1         0.0009      0.000892    0.0036      0.003573    0.01     0.03     |
|   compute_max_error_bound()           21        0.0083      0.000394    1.3041      0.062102    0.07     10.31    |
|   compute_output_dual_innerprods()    1         0.0048      0.004823    0.0198      0.019828    0.04     0.16     |
|   train_reduced_basis()               1         0.0108      0.010808    11.9772     11.977156   0.09     94.70    |
|   update_RB_system_matrices()         20        0.0591      0.002956    0.0591      0.002956    0.47     0.47     |
|   update_residual_terms()             20        0.1171      0.005856    0.2628      0.013140    0.93     2.08     |
|                                                                                                                   |
| RBEvaluation                                                                                                      |
|   clear()                             2         0.0006      0.000280    0.0006      0.000280    0.00     0.00     |
|   resize_data_structures()            1         0.0005      0.000549    0.0005      0.000549    0.00     0.00     |
|   write_offline_data_to_files()       1         0.0020      0.001981    0.0020      0.001981    0.02     0.02     |
|   write_out_basis_functions()         1         0.0000      0.000045    0.0706      0.070565    0.00     0.56     |
|   write_out_vectors()                 1         0.0396      0.039566    0.0705      0.070520    0.31     0.56     |
|                                                                                                                   |
| TransientRBConstruction                                                                                           |
|   enrich_RB_space()                   20        0.6043      0.030214    0.6043      0.030214    4.78     4.78     |
|   mass_matrix_scaled_matvec()         2000      0.2437      0.000122    0.2437      0.000122    1.93     1.93     |
|   set_error_temporal_data()           2020      0.5885      0.000291    0.5885      0.000291    4.65     4.65     |
|   truth_assembly()                    2000      6.7576      0.003379    7.0046      0.003502    53.43    55.38    |
|   truth_solve()                       20        0.5010      0.025050    9.5572      0.477860    3.96     75.57    |
|   update_RB_initial_condition_all_N() 20        0.0131      0.000655    0.0131      0.000655    0.10     0.10     |
|   update_RB_system_matrices()         20        0.0182      0.000911    0.0773      0.003867    0.14     0.61     |
|   update_residual_terms()             20        0.0780      0.003898    0.3860      0.019302    0.62     3.05     |
|                                                                                                                   |
| TransientRBEvaluation                                                                                             |
|   cache_online_residual_terms()       189       0.0157      0.000083    0.0157      0.000083    0.12     0.12     |
|   compute_residual_dual_norm()        18900     0.5972      0.000032    0.5972      0.000032    4.72     4.72     |
|   rb_solve()                          189       0.4952      0.002620    1.1126      0.005887    3.92     8.80     |
|   resize_data_structures()            1         0.0002      0.000235    0.0008      0.000785    0.00     0.01     |
|   write_offline_data_to_files()       1         0.0007      0.000733    0.0027      0.002716    0.01     0.02     |
 -------------------------------------------------------------------------------------------------------------------
| Totals:                               36791     12.6475                                         100.00            |
 -------------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example reduced_basis_ex3:
*  mpirun -np 12 example-devel -online_mode 0 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
***************************************************************
***************************************************************
* Running Example reduced_basis_ex3:
*  mpirun -np 12 example-devel -online_mode 1 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
***************************************************************
 
*** Warning, This code is untested, experimental, or likely to see future API changes: src/reduced_basis/rb_parametrized.C, line 41, compiled Jan 31 2013 at 21:51:32 ***
 EquationSystems
  n_systems()=1
   System #0, "RBConvectionDiffusion"
    Type "TransientRBConstruction"
    Variables="u" 
    Finite Element Types="LAGRANGE", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" 
    Approximation Orders="FIRST", "THIRD" 
    n_dofs()=676
    n_local_dofs()=68
    n_constrained_dofs()=100
    n_local_constrained_dofs()=16
    n_vectors()=3
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 7.99852
      Average Off-Processor Bandwidth <= 1.20118
      Maximum  On-Processor Bandwidth <= 11
      Maximum Off-Processor Bandwidth <= 9
    DofMap Constraints
      Number of DoF Constraints = 100
      Average DoF Constraint Length= 0
      Number of Node Constraints = 0

 Mesh Information:
  mesh_dimension()=2
  spatial_dimension()=3
  n_nodes()=676
    n_local_nodes()=68
  n_elem()=625
    n_local_elem()=51
    n_active_elem()=625
  n_subdomains()=1
  n_partitions()=12
  n_processors()=12
  n_threads()=1
  processor_id()=0

x_vel: 1.000000e+00
y_vel: 1.000000e+00

Error bound (absolute) at the final time is 0.0167619

************************************************************************************************************************
***             WIDEN YOUR WINDOW TO 120 CHARACTERS.  Use 'enscript -r -fCourier9' to print this document            ***
************************************************************************************************************************

---------------------------------------------- PETSc Performance Summary: ----------------------------------------------

/workspace/libmesh/examples/reduced_basis/reduced_basis_ex3/.libs/lt-example-devel on a intel-12. named hbar.ices.utexas.edu with 12 processors, by benkirk Thu Jan 31 22:17:35 2013
Using Petsc Release Version 3.3.0, Patch 2, Fri Jul 13 15:42:00 CDT 2012 

                         Max       Max/Min        Avg      Total 
Time (sec):           4.583e-01      1.00062   4.583e-01
Objects:              5.700e+01      1.00000   5.700e+01
Flops:                2.720e+03      1.78947   2.253e+03  2.704e+04
Flops/sec:            5.935e+03      1.78947   4.917e+03  5.900e+04
MPI Messages:         7.800e+01      3.25000   5.000e+01  6.000e+02
MPI Message Lengths:  2.370e+03      2.38431   3.191e+01  1.914e+04
MPI Reductions:       1.390e+02      1.00000

Flop counting convention: 1 flop = 1 real number operation of type (multiply/divide/add/subtract)
                            e.g., VecAXPY() for real vectors of length N --> 2N flops
                            and VecAXPY() for complex vectors of length N --> 8N flops

Summary of Stages:   ----- Time ------  ----- Flops -----  --- Messages ---  -- Message Lengths --  -- Reductions --
                        Avg     %Total     Avg     %Total   counts   %Total     Avg         %Total   counts   %Total 
 0:      Main Stage: 4.5824e-01 100.0%  2.7040e+04 100.0%  6.000e+02 100.0%  3.191e+01      100.0%  1.380e+02  99.3% 

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

VecCopy                3 1.0 1.8120e-05 1.9 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecSet                29 1.0 1.9550e-05 1.8 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecAXPY               20 1.0 5.1260e-05 1.8 2.72e+03 1.8 0.0e+00 0.0e+00 0.0e+00  0100  0  0  0   0100  0  0  0   528
VecAssemblyBegin      23 1.0 7.1752e-0315.6 0.00e+00 0.0 0.0e+00 0.0e+00 6.9e+01  1  0  0  0 50   1  0  0  0 50     0
VecAssemblyEnd        23 1.0 4.5061e-05 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecScatterBegin        2 1.0 9.2983e-05 1.1 0.00e+00 0.0 1.0e+02 8.1e+01 0.0e+00  0  0 17 42  0   0  0 17 42  0     0
VecScatterEnd          2 1.0 2.5988e-05 1.8 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatZeroEntries         2 1.0 1.8120e-05 1.8 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
------------------------------------------------------------------------------------------------------------------------

Memory usage is given in bytes:

Object Type          Creations   Destructions     Memory  Descendants' Mem.
Reports information only for process 0.

--- Event Stage 0: Main Stage

              Vector    33             33        65080     0
      Vector Scatter     5              5         5180     0
           Index Set    10             10         7800     0
   IS L to G Mapping     5              5         2820     0
              Matrix     3              3        16592     0
              Viewer     1              0            0     0
========================================================================================================================
Average time to get PetscTime(): 0
Average time for MPI_Barrier(): 4.76837e-06
Average time for zero size MPI_Send(): 1.34905e-05
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
| Time:           Thu Jan 31 22:17:35 2013                                                                             |
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
 --------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=0.505343, Active time=0.436629                                               |
 --------------------------------------------------------------------------------------------------------------
| Event                            nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                            w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|--------------------------------------------------------------------------------------------------------------|
|                                                                                                              |
|                                                                                                              |
| DofMap                                                                                                       |
|   add_neighbors_to_send_list()   1         0.0031      0.003147    0.0051      0.005051    0.72     1.16     |
|   build_sparsity()               1         0.0021      0.002109    0.0066      0.006616    0.48     1.52     |
|   create_dof_constraints()       1         0.0117      0.011748    0.0463      0.046269    2.69     10.60    |
|   distribute_dofs()              1         0.0111      0.011062    0.0367      0.036725    2.53     8.41     |
|   dof_indices()                  813       0.0445      0.000055    0.0445      0.000055    10.20    10.20    |
|   prepare_send_list()            1         0.0000      0.000037    0.0000      0.000037    0.01     0.01     |
|   reinit()                       1         0.0236      0.023610    0.0236      0.023610    5.41     5.41     |
|                                                                                                              |
| EquationSystems                                                                                              |
|   build_solution_vector()        2         0.0017      0.000866    0.0084      0.004184    0.40     1.92     |
|                                                                                                              |
| ExodusII_IO                                                                                                  |
|   write_nodal_data()             2         0.0126      0.006306    0.0126      0.006306    2.89     2.89     |
|                                                                                                              |
| Mesh                                                                                                         |
|   find_neighbors()               1         0.0139      0.013879    0.0144      0.014353    3.18     3.29     |
|   renumber_nodes_and_elem()      2         0.0006      0.000290    0.0006      0.000290    0.13     0.13     |
|                                                                                                              |
| MeshCommunication                                                                                            |
|   assign_global_indices()        1         0.0261      0.026107    0.0838      0.083813    5.98     19.20    |
|   compute_hilbert_indices()      2         0.0108      0.005409    0.0108      0.005409    2.48     2.48     |
|   find_global_indices()          2         0.0044      0.002193    0.0199      0.009953    1.00     4.56     |
|   parallel_sort()                2         0.0029      0.001475    0.0035      0.001734    0.68     0.79     |
|                                                                                                              |
| MeshOutput                                                                                                   |
|   write_equation_systems()       2         0.0002      0.000085    0.0214      0.010695    0.04     4.90     |
|                                                                                                              |
| MeshTools::Generation                                                                                        |
|   build_cube()                   1         0.0032      0.003178    0.0032      0.003178    0.73     0.73     |
|                                                                                                              |
| MetisPartitioner                                                                                             |
|   partition()                    1         0.0495      0.049543    0.0586      0.058633    11.35    13.43    |
|                                                                                                              |
| Parallel                                                                                                     |
|   allgather()                    14        0.0274      0.001960    0.0284      0.002027    6.28     6.50     |
|   barrier()                      1         0.0710      0.070987    0.0710      0.070987    16.26    16.26    |
|   broadcast()                    13        0.0002      0.000018    0.0002      0.000015    0.05     0.05     |
|   max(bool)                      1         0.0000      0.000007    0.0000      0.000007    0.00     0.00     |
|   max(scalar)                    135       0.0038      0.000028    0.0038      0.000028    0.86     0.86     |
|   max(vector)                    30        0.0009      0.000030    0.0031      0.000102    0.20     0.70     |
|   min(bool)                      156       0.0035      0.000023    0.0035      0.000023    0.81     0.81     |
|   min(scalar)                    129       0.0359      0.000279    0.0359      0.000279    8.23     8.23     |
|   min(vector)                    30        0.0010      0.000032    0.0049      0.000163    0.22     1.12     |
|   probe()                        224       0.0041      0.000018    0.0041      0.000018    0.94     0.94     |
|   receive()                      202       0.0012      0.000006    0.0053      0.000026    0.27     1.20     |
|   send()                         202       0.0006      0.000003    0.0006      0.000003    0.14     0.14     |
|   send_receive()                 184       0.0014      0.000007    0.0073      0.000040    0.31     1.67     |
|   sum()                          31        0.0100      0.000321    0.0357      0.001151    2.28     8.17     |
|                                                                                                              |
| Parallel::Request                                                                                            |
|   wait()                         202       0.0003      0.000002    0.0003      0.000002    0.08     0.08     |
|                                                                                                              |
| Partitioner                                                                                                  |
|   set_node_processor_ids()       1         0.0016      0.001588    0.0024      0.002359    0.36     0.54     |
|   set_parent_processor_ids()     1         0.0012      0.001162    0.0012      0.001162    0.27     0.27     |
|                                                                                                              |
| RBConstruction                                                                                               |
|   clear()                        3         0.0005      0.000168    0.0005      0.000168    0.12     0.12     |
|   load_basis_function()          1         0.0001      0.000096    0.0001      0.000096    0.02     0.02     |
|                                                                                                              |
| RBEvaluation                                                                                                 |
|   clear()                        2         0.0001      0.000045    0.0001      0.000045    0.02     0.02     |
|   read_in_basis_functions()      1         0.0001      0.000051    0.1957      0.195667    0.01     44.81    |
|   read_in_vectors()              1         0.0218      0.021796    0.1956      0.195616    4.99     44.80    |
|   read_offline_data_from_files() 1         0.0018      0.001839    0.0028      0.002827    0.42     0.65     |
|   resize_data_structures()       1         0.0007      0.000692    0.0007      0.000692    0.16     0.16     |
|                                                                                                              |
| TransientRBConstruction                                                                                      |
|   load_rb_solution()             1         0.0004      0.000385    0.0004      0.000385    0.09     0.09     |
|                                                                                                              |
| TransientRBEvaluation                                                                                        |
|   cache_online_residual_terms()  1         0.0002      0.000230    0.0002      0.000230    0.05     0.05     |
|   compute_residual_dual_norm()   100       0.0047      0.000047    0.0047      0.000047    1.08     1.08     |
|   rb_solve()                     1         0.0191      0.019113    0.0241      0.024060    4.38     5.51     |
|   read_offline_data_from_files() 1         0.0006      0.000579    0.0034      0.003406    0.13     0.78     |
|   resize_data_structures()       1         0.0003      0.000296    0.0010      0.000988    0.07     0.23     |
 --------------------------------------------------------------------------------------------------------------
| Totals:                          2508      0.4366                                          100.00            |
 --------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example reduced_basis_ex3:
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
