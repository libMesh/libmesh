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
<br><br><br> <h1> The source file assembly.h with comments: </h1> 
<div class = "comment">
</div>

<div class ="fragment">
<pre>
        #ifndef __assembly_h__
        #define __assembly_h__
        
        #if defined(LIBMESH_HAVE_SLEPC) && defined(LIBMESH_HAVE_GLPK)
        
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
        using libMesh::RBTheta;
        using libMesh::Real;
        using libMesh::RealGradient;
        
</pre>
</div>
<div class = "comment">
Functors for the parameter-dependent part of the affine decomposition of the PDE
The RHS and outputs just require a constant value of 1, so use a default RBTheta object there
</div>

<div class ="fragment">
<pre>
        struct ThetaA0 : RBTheta { virtual Number evaluate(const RBParameters& mu) { return mu.get_value("mu_0"); } };
        struct ThetaA1 : RBTheta { virtual Number evaluate(const RBParameters& mu) { return mu.get_value("mu_1"); } };
        struct ThetaA2 : RBTheta { virtual Number evaluate(const RBParameters& mu) { return mu.get_value("mu_2"); } };
        
        struct B : ElemAssembly
        {
</pre>
</div>
<div class = "comment">
Assemble the H1 inner product. This will be used as the inner product
for this problem.
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
                  c.get_elem_jacobian()(i,j) += JxW[qp] * (phi[j][qp]*phi[i][qp] + dphi[j][qp]*dphi[i][qp]);
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
        
            Real min_x=0.,max_x=0.5;
        
            Point centroid = c.elem-&gt;centroid();
            if( (min_x &lt;= centroid(0)) && (centroid(0) &lt; max_x) )
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
        
            Real min_x=0.5, max_x=1.;
        
            Point centroid = c.elem-&gt;centroid();
            if( (min_x &lt;= centroid(0)) && (centroid(0) &lt;= max_x) )
            for (unsigned int qp=0; qp != n_qpoints; qp++)
              for (unsigned int i=0; i != n_u_dofs; i++)
                for (unsigned int j=0; j != n_u_dofs; j++)
                  c.get_elem_jacobian()(i,j) += JxW[qp] * dphi[j][qp]*dphi[i][qp];
          }
        };
        
        struct A2 : ElemAssembly
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
                  c.get_elem_jacobian()(i,j) -= JxW[qp] *dphi[i][qp](0)*phi[j][qp];
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
        struct Ex02RBThetaExpansion : RBThetaExpansion
        {
        
          /**
           * Constructor.
           */
          Ex02RBThetaExpansion()
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
        struct Ex02RBAssemblyExpansion : RBAssemblyExpansion
        {
        
          /**
           * Constructor.
           */
          Ex02RBAssemblyExpansion()
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
          B B_assembly;
          A0 A0_assembly;
          A1 A1_assembly;
          A2 A2_assembly;
          F0 F0_assembly;
          OutputAssembly L0;
          OutputAssembly L1;
          OutputAssembly L2;
          OutputAssembly L3;
        };
        
        #endif // LIBMESH_HAVE_SLEPC && LIBMESH_HAVE_GLPK
        
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
        
        #if defined(LIBMESH_HAVE_SLEPC) && defined(LIBMESH_HAVE_GLPK)
        
        #include "libmesh/rb_construction.h"
        #include "libmesh/rb_scm_construction.h"
        #include "libmesh/fe_base.h"
        
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
        using libMesh::RBSCMConstruction;
        using libMesh::RBEvaluation;
        using libMesh::RBSCMEvaluation;
        using libMesh::RBParameters;
        using libMesh::RBThetaExpansion;
        using libMesh::RBAssemblyExpansion;
        using libMesh::AutoPtr;
        using libMesh::DirichletBoundary;
        using libMesh::Real;
        
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
A simple subclass of RBEvaluation. We also store the theta expansion object
for the affine expansion of the PDE as a member variable.
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
            set_rb_theta_expansion(ex02_rb_theta_expansion);
          }
        
          /**
           * We override get_stability_lower_bound so that it calls rb_scm_eval to return
           * a parameter-dependent lower bound for the coercivity constant.
           */
          virtual Real get_stability_lower_bound()
          { 
            rb_scm_eval-&gt;set_parameters( get_parameters() );
            return rb_scm_eval-&gt;get_SCM_LB() ;
          }
        
          /**
           * Pointer to the SCM object that will provide our coercivity constant lower bound.
           */
          RBSCMEvaluation* rb_scm_eval;
          
          /**
           * The object that stores the "theta" expansion of the parameter dependent PDE,
           * i.e. the set of parameter-dependent functions in the affine expansion of the PDE.
           */
          Ex02RBThetaExpansion ex02_rb_theta_expansion;
        
        };
        
</pre>
</div>
<div class = "comment">
A simple subclass of RBConstruction, which initializes libMesh-related data such
as the number of variables and their finite element type. We also store the objects
that define the affine expansion of the PDE as member variables.
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
          virtual ~SimpleRBConstruction () {}
        
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
            set_rb_assembly_expansion(ex02_rb_assembly_expansion);
        
</pre>
</div>
<div class = "comment">
We need to define an inner product matrix for this problem
</div>

<div class ="fragment">
<pre>
            set_inner_product_assembly(ex02_rb_assembly_expansion.B_assembly);
          }
        
          /**
           * Pre-request all relevant element data. (This is not essential, but it
           * allows libMesh to cache data and hence run faster.)
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
          Ex02RBAssemblyExpansion ex02_rb_assembly_expansion;
        
          /**
           * The object that defines which degrees of freedom are on a Dirichlet boundary.
           */
          AutoPtr&lt;DirichletBoundary&gt; dirichlet_bc;
        
        };
        
        #endif // LIBMESH_HAVE_SLEPC && LIBMESH_HAVE_GLPK
        
        #endif
</pre>
</div>

<a name="comments"></a> 
<br><br><br> <h1> The source file reduced_basis_ex2.C with comments: </h1> 
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
        #include &lt;cstdlib&gt; // *must* precede &lt;cmath&gt; for proper std:abs() on PGI, Sun Studio CC
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
<br><br><br> <h1> The source file assembly.h without comments: </h1> 
<pre> 
  #ifndef __assembly_h__
  #define __assembly_h__
  
  #<B><FONT COLOR="#A020F0">if</FONT></B> defined(LIBMESH_HAVE_SLEPC) &amp;&amp; defined(LIBMESH_HAVE_GLPK)
  
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
  
  using libMesh::ElemAssembly;
  using libMesh::FEInterface;
  using libMesh::FEMContext;
  using libMesh::Number;
  using libMesh::Point;
  using libMesh::RBTheta;
  using libMesh::Real;
  using libMesh::RealGradient;
  
  <B><FONT COLOR="#228B22">struct</FONT></B> ThetaA0 : RBTheta { <B><FONT COLOR="#228B22">virtual</FONT></B> Number evaluate(<B><FONT COLOR="#228B22">const</FONT></B> RBParameters&amp; mu) { <B><FONT COLOR="#A020F0">return</FONT></B> mu.get_value(<B><FONT COLOR="#BC8F8F">&quot;mu_0&quot;</FONT></B>); } };
  <B><FONT COLOR="#228B22">struct</FONT></B> ThetaA1 : RBTheta { <B><FONT COLOR="#228B22">virtual</FONT></B> Number evaluate(<B><FONT COLOR="#228B22">const</FONT></B> RBParameters&amp; mu) { <B><FONT COLOR="#A020F0">return</FONT></B> mu.get_value(<B><FONT COLOR="#BC8F8F">&quot;mu_1&quot;</FONT></B>); } };
  <B><FONT COLOR="#228B22">struct</FONT></B> ThetaA2 : RBTheta { <B><FONT COLOR="#228B22">virtual</FONT></B> Number evaluate(<B><FONT COLOR="#228B22">const</FONT></B> RBParameters&amp; mu) { <B><FONT COLOR="#A020F0">return</FONT></B> mu.get_value(<B><FONT COLOR="#BC8F8F">&quot;mu_2&quot;</FONT></B>); } };
  
  <B><FONT COLOR="#228B22">struct</FONT></B> B : ElemAssembly
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
            c.get_elem_jacobian()(i,j) += JxW[qp] * (phi[j][qp]*phi[i][qp] + dphi[j][qp]*dphi[i][qp]);
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
  
      Real min_x=0.,max_x=0.5;
  
      Point centroid = c.elem-&gt;centroid();
      <B><FONT COLOR="#A020F0">if</FONT></B>( (min_x &lt;= centroid(0)) &amp;&amp; (centroid(0) &lt; max_x) )
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
  
      <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;RealGradient&gt; &gt;&amp; dphi =
        c.element_fe_var[u_var]-&gt;get_dphi();
  
      <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_u_dofs = c.dof_indices_var[u_var].size();
  
      <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_qpoints = (c.get_element_qrule())-&gt;n_points();
  
      Real min_x=0.5, max_x=1.;
  
      Point centroid = c.elem-&gt;centroid();
      <B><FONT COLOR="#A020F0">if</FONT></B>( (min_x &lt;= centroid(0)) &amp;&amp; (centroid(0) &lt;= max_x) )
      <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> qp=0; qp != n_qpoints; qp++)
        <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i != n_u_dofs; i++)
          <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j=0; j != n_u_dofs; j++)
            c.get_elem_jacobian()(i,j) += JxW[qp] * dphi[j][qp]*dphi[i][qp];
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
            c.get_elem_jacobian()(i,j) -= JxW[qp] *dphi[i][qp](0)*phi[j][qp];
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
  
  <B><FONT COLOR="#228B22">struct</FONT></B> Ex02RBThetaExpansion : RBThetaExpansion
  {
  
    <I><FONT COLOR="#B22222">/**
     * Constructor.
     */</FONT></I>
    Ex02RBThetaExpansion()
    {
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
  
  <B><FONT COLOR="#228B22">struct</FONT></B> Ex02RBAssemblyExpansion : RBAssemblyExpansion
  {
  
    <I><FONT COLOR="#B22222">/**
     * Constructor.
     */</FONT></I>
    Ex02RBAssemblyExpansion()
      :
      L0(0.72,0.88,0.72,0.88), <I><FONT COLOR="#B22222">// We make sure these output regions conform to the mesh
</FONT></I>      L1(0.12,0.28,0.72,0.88),
      L2(0.12,0.28,0.12,0.28),
      L3(0.72,0.88,0.12,0.28)
    {
      attach_A_assembly(&amp;A0_assembly); <I><FONT COLOR="#B22222">// Attach the lhs assembly
</FONT></I>      attach_A_assembly(&amp;A1_assembly);
      attach_A_assembly(&amp;A2_assembly);
      
      attach_F_assembly(&amp;F0_assembly); <I><FONT COLOR="#B22222">// Attach the rhs assembly
</FONT></I>      
      attach_output_assembly(&amp;L0);       <I><FONT COLOR="#B22222">// Attach output 0 assembly
</FONT></I>      attach_output_assembly(&amp;L1);       <I><FONT COLOR="#B22222">// Attach output 1 assembly
</FONT></I>      attach_output_assembly(&amp;L2);       <I><FONT COLOR="#B22222">// Attach output 2 assembly
</FONT></I>      attach_output_assembly(&amp;L3);       <I><FONT COLOR="#B22222">// Attach output 3 assembly
</FONT></I>    }
  
    B B_assembly;
    A0 A0_assembly;
    A1 A1_assembly;
    A2 A2_assembly;
    F0 F0_assembly;
    OutputAssembly L0;
    OutputAssembly L1;
    OutputAssembly L2;
    OutputAssembly L3;
  };
  
  #endif <I><FONT COLOR="#B22222">// LIBMESH_HAVE_SLEPC &amp;&amp; LIBMESH_HAVE_GLPK
</FONT></I>  
  #endif
  
</pre> 
<a name="nocomments"></a> 
<br><br><br> <h1> The source file rb_classes.h without comments: </h1> 
<pre> 
    
    
  
  #ifndef __rb_classes_h__
  #define __rb_classes_h__
  
  #<B><FONT COLOR="#A020F0">if</FONT></B> defined(LIBMESH_HAVE_SLEPC) &amp;&amp; defined(LIBMESH_HAVE_GLPK)
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/rb_construction.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/rb_scm_construction.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/fe_base.h&quot;</FONT></B>
  
  using libMesh::EquationSystems;
  using libMesh::FEMContext;
  using libMesh::RBConstruction;
  using libMesh::RBSCMConstruction;
  using libMesh::RBEvaluation;
  using libMesh::RBSCMEvaluation;
  using libMesh::RBParameters;
  using libMesh::RBThetaExpansion;
  using libMesh::RBAssemblyExpansion;
  using libMesh::AutoPtr;
  using libMesh::DirichletBoundary;
  using libMesh::Real;
  
  #include <B><FONT COLOR="#BC8F8F">&quot;assembly.h&quot;</FONT></B>
  
  <B><FONT COLOR="#228B22">class</FONT></B> SimpleRBEvaluation : <B><FONT COLOR="#228B22">public</FONT></B> RBEvaluation
  {
  <B><FONT COLOR="#228B22">public</FONT></B>:
  
    <I><FONT COLOR="#B22222">/**
     * Constructor. Just set the theta expansion.
     */</FONT></I>
    SimpleRBEvaluation()
    {
      set_rb_theta_expansion(ex02_rb_theta_expansion);
    }
  
    <I><FONT COLOR="#B22222">/**
     * We override get_stability_lower_bound so that it calls rb_scm_eval to return
     * a parameter-dependent lower bound for the coercivity constant.
     */</FONT></I>
    <B><FONT COLOR="#228B22">virtual</FONT></B> Real get_stability_lower_bound()
    { 
      rb_scm_eval-&gt;set_parameters( get_parameters() );
      <B><FONT COLOR="#A020F0">return</FONT></B> rb_scm_eval-&gt;get_SCM_LB() ;
    }
  
    <I><FONT COLOR="#B22222">/**
     * Pointer to the SCM object that will provide our coercivity constant lower bound.
     */</FONT></I>
    RBSCMEvaluation* rb_scm_eval;
    
    <I><FONT COLOR="#B22222">/**
     * The object that stores the &quot;theta&quot; expansion of the parameter dependent PDE,
     * i.e. the set of parameter-dependent functions in the affine expansion of the PDE.
     */</FONT></I>
    Ex02RBThetaExpansion ex02_rb_theta_expansion;
  
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
    <B><FONT COLOR="#228B22">virtual</FONT></B> ~SimpleRBConstruction () {}
  
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
      u_var = <B><FONT COLOR="#A020F0">this</FONT></B>-&gt;add_variable (<B><FONT COLOR="#BC8F8F">&quot;u&quot;</FONT></B>, FIRST);
  
      dirichlet_bc = build_zero_dirichlet_boundary_object();
      
      dirichlet_bc-&gt;b.insert(3);
      dirichlet_bc-&gt;variables.push_back(u_var);
      
      get_dof_map().add_dirichlet_boundary(*dirichlet_bc);
  
      <B><FONT COLOR="#5F9EA0">Parent</FONT></B>::init_data();
      
      set_rb_assembly_expansion(ex02_rb_assembly_expansion);
  
      set_inner_product_assembly(ex02_rb_assembly_expansion.B_assembly);
    }
  
    <I><FONT COLOR="#B22222">/**
     * Pre-request all relevant element data. (This is not essential, but it
     * allows libMesh to cache data and hence run faster.)
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
    Ex02RBAssemblyExpansion ex02_rb_assembly_expansion;
  
    <I><FONT COLOR="#B22222">/**
     * The object that defines which degrees of freedom are on a Dirichlet boundary.
     */</FONT></I>
    AutoPtr&lt;DirichletBoundary&gt; dirichlet_bc;
  
  };
  
  #endif <I><FONT COLOR="#B22222">// LIBMESH_HAVE_SLEPC &amp;&amp; LIBMESH_HAVE_GLPK
</FONT></I>  
  #endif
</pre> 
<a name="nocomments"></a> 
<br><br><br> <h1> The source file reduced_basis_ex2.C without comments: </h1> 
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
  #include &lt;cstdlib&gt; <I><FONT COLOR="#B22222">// *must* precede &lt;cmath&gt; for proper std:abs() on PGI, Sun Studio CC
</FONT></I>  #include &lt;cmath&gt;
  #include &lt;set&gt;
  
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
***************************************************************
* Running Example reduced_basis_ex2:
*  mpirun -np 12 example-devel -online_mode 0 -eps_type lapack -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
***************************************************************
 
*** Warning, This code is untested, experimental, or likely to see future API changes: src/reduced_basis/rb_parametrized.C, line 41, compiled Jan 31 2013 at 21:51:32 ***
 EquationSystems
  n_systems()=2
   System #0, "RBConvectionDiffusion"
    Type "RBConstruction"
    Variables="u" 
    Finite Element Types="LAGRANGE", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" 
    Approximation Orders="FIRST", "THIRD" 
    n_dofs()=676
    n_local_dofs()=68
    n_constrained_dofs()=26
    n_local_constrained_dofs()=6
    n_vectors()=1
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 7.99852
      Average Off-Processor Bandwidth <= 1.20118
      Maximum  On-Processor Bandwidth <= 11
      Maximum Off-Processor Bandwidth <= 9
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
    n_local_dofs()=68
    n_constrained_dofs()=0
    n_local_constrained_dofs()=0
    n_vectors()=0
    n_matrices()=2
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 7.99852
      Average Off-Processor Bandwidth <= 1.20118
      Maximum  On-Processor Bandwidth <= 11
      Maximum Off-Processor Bandwidth <= 9
    DofMap Constraints
      Number of DoF Constraints = 0
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


B_min(0) = -1.19047e-15
B_max(0) = 0.999932

B_min(1) = -9.86284e-16
B_max(1) = 0.999933

B_min(2) = -0.380786
B_max(2) = 5.42317e-17

SCM: Added mu = (0.69213,0.358072,0.0138549)

Stability constant for C_J(0) = 0.312012

SCM iteration 0, max_SCM_error = 0.956432

SCM: Added mu = (0.152571,0.950732,0.0958166)

-----------------------------------


Stability constant for C_J(1) = 0.0962708

SCM iteration 1, max_SCM_error = 0.521896

SCM: Added mu = (0.130352,0.12171,0.0598661)

-----------------------------------


Stability constant for C_J(2) = 0.0710684

SCM iteration 2, max_SCM_error = 0.275696

SCM: Added mu = (0.732283,0.102329,0.0756287)

-----------------------------------


Stability constant for C_J(3) = 0.0742499

SCM iteration 3, max_SCM_error = 0.261343

SCM: Added mu = (0.917311,0.961874,0.011123)

-----------------------------------


Stability constant for C_J(4) = 0.655301

SCM iteration 4, max_SCM_error = 0.202518

SCM: Added mu = (0.116516,0.210237,0.0775824)

-----------------------------------


Stability constant for C_J(5) = 0.0683307

SCM iteration 5, max_SCM_error = 0.125482

SCM: Added mu = (0.159746,0.11326,0.0485068)

-----------------------------------


Stability constant for C_J(6) = 0.0800011

SCM iteration 6, max_SCM_error = 0.111669

SCM: Added mu = (0.320456,0.115153,0.0760262)

-----------------------------------


Stability constant for C_J(7) = 0.0827318

SCM iteration 7, max_SCM_error = 0.0983686

SCM tolerance of 0.1 reached.

Compute output dual inner products
output_dual_innerprods[0][0] = 0.839698
output_dual_innerprods[1][0] = 0.318298
output_dual_innerprods[2][0] = 0.318298
output_dual_innerprods[3][0] = 0.839698

---- Performing Greedy basis enrichment ----

---- Basis dimension: 0 ----
Performing RB solves on training set
Maximum (absolute) error bound is 7.14651

Performing truth solve at parameter:
mu_0: 1.106337e-01
mu_1: 1.858370e-01
mu_2: 5.639861e-02

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 1 ----
Performing RB solves on training set
Maximum (absolute) error bound is 8.27415

Performing truth solve at parameter:
mu_0: 1.112930e-01
mu_1: 7.834330e-01
mu_2: 1.870024e-02

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 2 ----
Performing RB solves on training set
Maximum (absolute) error bound is 4.13511

Performing truth solve at parameter:
mu_0: 7.322826e-01
mu_1: 1.023286e-01
mu_2: 7.562866e-02

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 3 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.475602

Performing truth solve at parameter:
mu_0: 7.665338e-01
mu_1: 1.186039e-01
mu_2: 1.243793e-02

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 4 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.143949

Performing truth solve at parameter:
mu_0: 1.525709e-01
mu_1: 9.507324e-01
mu_2: 9.581664e-02

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 5 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.00803168

Performing truth solve at parameter:
mu_0: 9.834567e-01
mu_1: 1.036751e-01
mu_2: 2.605954e-02

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 6 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.00354437

Performing truth solve at parameter:
mu_0: 1.575982e-01
mu_1: 8.958385e-01
mu_2: 6.370547e-02

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 7 ----
Performing RB solves on training set
Maximum (absolute) error bound is 4.70014e-05

Performing truth solve at parameter:
mu_0: 3.204558e-01
mu_1: 1.151532e-01
mu_2: 7.602619e-02

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 8 ----
Performing RB solves on training set
Maximum (absolute) error bound is 7.90815e-06

Specified error tolerance reached.
************************************************************************************************************************
***             WIDEN YOUR WINDOW TO 120 CHARACTERS.  Use 'enscript -r -fCourier9' to print this document            ***
************************************************************************************************************************

---------------------------------------------- PETSc Performance Summary: ----------------------------------------------

/workspace/libmesh/examples/reduced_basis/reduced_basis_ex2/.libs/lt-example-devel on a intel-12. named hbar.ices.utexas.edu with 12 processors, by benkirk Thu Jan 31 22:16:52 2013
Using Petsc Release Version 3.3.0, Patch 2, Fri Jul 13 15:42:00 CDT 2012 

                         Max       Max/Min        Avg      Total 
Time (sec):           8.002e+00      1.00000   8.002e+00
Objects:              1.024e+04      1.00000   1.024e+04
Flops:                2.788e+07      1.86785   2.234e+07  2.681e+08
Flops/sec:            3.484e+06      1.86785   2.792e+06  3.351e+07
MPI Messages:         2.719e+04      3.00304   1.812e+04  2.174e+05
MPI Message Lengths:  1.516e+06      2.27314   5.676e+01  1.234e+07
MPI Reductions:       3.060e+04      1.00000

Flop counting convention: 1 flop = 1 real number operation of type (multiply/divide/add/subtract)
                            e.g., VecAXPY() for real vectors of length N --> 2N flops
                            and VecAXPY() for complex vectors of length N --> 8N flops

Summary of Stages:   ----- Time ------  ----- Flops -----  --- Messages ---  -- Message Lengths --  -- Reductions --
                        Avg     %Total     Avg     %Total   counts   %Total     Avg         %Total   counts   %Total 
 0:      Main Stage: 8.0021e+00 100.0%  2.6813e+08 100.0%  2.174e+05 100.0%  5.676e+01      100.0%  3.060e+04 100.0% 

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

STSetUp               14 1.0 7.7524e-02 9.6 0.00e+00 0.0 0.0e+00 0.0e+00 2.8e+01  0  0  0  0  0   0  0  0  0  0     0
EPSSetUp              14 1.0 4.9376e-01 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 1.9e+04  6  0  0  0 61   6  0  0  0 61     0
EPSSolve              14 1.0 5.4305e+00 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00 66  0  0  0  0  66  0  0  0  0     0
DSSolve               14 1.0 5.3821e+00 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00 65  0  0  0  0  65  0  0  0  0     0
DSVectors             14 1.0 2.8543e-02 2.5 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
DSOther               14 1.0 2.1785e-02 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecDot               753 1.0 2.4075e-02 3.9 1.02e+05 1.8 0.0e+00 0.0e+00 7.5e+02  0  0  0  0  2   0  0  0  0  2    42
VecMDot             3426 1.0 5.7369e-02 1.3 7.02e+06 1.8 0.0e+00 0.0e+00 3.4e+03  1 26  0  0 11   1 26  0  0 11  1215
VecNorm             3596 1.0 3.6439e-02 1.1 4.89e+05 1.8 0.0e+00 0.0e+00 3.6e+03  0  2  0  0 12   0  2  0  0 12   133
VecScale            3599 1.0 1.8098e-03 1.5 2.44e+05 1.8 0.0e+00 0.0e+00 0.0e+00  0  1  0  0  0   0  1  0  0  0  1341
VecCopy              271 1.0 1.6165e-04 1.4 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecSet              4863 1.0 1.9207e-03 1.4 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecAXPY              307 1.0 2.8968e-04 1.5 4.18e+04 1.8 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0  1433
VecMAXPY            3559 1.0 5.5616e-03 1.5 7.54e+06 1.8 0.0e+00 0.0e+00 0.0e+00  0 28  0  0  0   0 28  0  0  0 13480
VecAssemblyBegin     247 1.0 3.7147e-02 4.0 0.00e+00 0.0 2.4e+02 8.1e+01 7.4e+02  0  0  0  0  2   0  0  0  0  2     0
VecAssemblyEnd       247 1.0 2.0051e-04 1.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecScatterBegin     4299 1.0 1.5754e-02 1.8 0.00e+00 0.0 2.1e+05 5.4e+01 0.0e+00  0  0 95 91  0   0  0 95 91  0     0
VecScatterEnd       4299 1.0 3.0859e-02 4.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecNormalize        3559 1.0 3.8416e-02 1.1 7.26e+05 1.8 0.0e+00 0.0e+00 3.6e+03  0  3  0  0 12   0  3  0  0 12   188
MatMult             3559 1.0 4.2127e-02 1.9 3.76e+06 1.6 1.7e+05 5.4e+01 0.0e+00  0 14 79 75  0   0 14 79 75  0   919
MatMultAdd           689 1.0 1.0989e-02 2.4 7.74e+05 1.6 3.3e+04 5.4e+01 0.0e+00  0  3 15 14  0   0  3 15 14  0   724
MatSolve            3596 1.0 1.1095e-02 2.4 7.75e+06 2.2 0.0e+00 0.0e+00 0.0e+00  0 25  0  0  0   0 25  0  0  0  6026
MatLUFactorNum        18 1.0 7.2145e-04 1.9 1.54e+05 2.6 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0  1670
MatILUFactorSym       18 1.0 1.4954e-03 1.6 0.00e+00 0.0 0.0e+00 0.0e+00 5.4e+01  0  0  0  0  0   0  0  0  0  0     0
MatConvert            28 1.0 7.4607e-02 1.7 0.00e+00 0.0 0.0e+00 0.0e+00 5.6e+01  1  0  0  0  0   1  0  0  0  0     0
MatAssemblyBegin    1068 1.0 4.1622e-01 3.0 0.00e+00 0.0 4.0e+03 2.5e+02 2.0e+03  3  0  2  8  6   3  0  2  8  6     0
MatAssemblyEnd      1068 1.0 5.2363e-02 3.7 0.00e+00 0.0 6.5e+03 1.5e+01 5.6e+02  0  0  3  1  2   0  0  3  1  2     0
MatGetRow          22960 1.1 2.4874e-03 1.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatGetRowIJ           18 1.0 1.6928e-05 3.5 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatGetSubMatrice      28 1.0 3.7783e-02 3.6 0.00e+00 0.0 0.0e+00 0.0e+00 2.2e+02  0  0  0  0  1   0  0  0  0  1     0
MatGetOrdering        18 1.0 5.6601e-04 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 3.6e+01  0  0  0  0  0   0  0  0  0  0     0
MatZeroEntries        72 1.0 1.2922e-04 1.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatAXPY               35 1.0 2.4285e-02 1.0 0.00e+00 0.0 3.4e+03 1.6e+01 5.6e+02  0  0  2  0  2   0  0  2  0  2     0
KSPGMRESOrthog      3426 1.0 6.5191e-02 1.3 1.41e+07 1.8 0.0e+00 0.0e+00 3.4e+03  1 52  0  0 11   1 52  0  0 11  2148
KSPSetUp              69 1.0 3.9148e-04 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPSolve              37 1.0 1.8056e-01 1.0 2.70e+07 1.9 1.7e+05 5.4e+01 7.1e+03  2 97 79 75 23   2 97 79 75 23  1435
PCSetUp               50 1.0 7.0183e-03 1.2 1.54e+05 2.6 0.0e+00 0.0e+00 9.2e+01  0  0  0  0  0   0  0  0  0  0   172
PCSetUpOnBlocks       37 1.0 5.7864e-03 1.2 1.54e+05 2.6 0.0e+00 0.0e+00 9.0e+01  0  0  0  0  0   0  0  0  0  0   208
PCApply             3596 1.0 4.2347e-02 1.2 7.75e+06 2.2 0.0e+00 0.0e+00 0.0e+00  0 25  0  0  0   0 25  0  0  0  1579
------------------------------------------------------------------------------------------------------------------------

Memory usage is given in bytes:

Object Type          Creations   Destructions     Memory  Descendants' Mem.
Reports information only for process 0.

--- Event Stage 0: Main Stage

           Container    14             14         7672     0
  Spectral Transform     1              1          736     0
 Eigenproblem Solver     1              1       328740     0
       Inner product     1              1          624     0
       Direct solver     1              1    284576440     0
              Vector  9461           9461     14263912     0
      Vector Scatter    72             72        74592     0
           Index Set   338            338       260832     0
   IS L to G Mapping     2              2         1128     0
              Matrix   340            340     98692760     0
         PetscRandom     1              1          608     0
       Krylov Solver     3              3        20432     0
      Preconditioner     3              3         2640     0
              Viewer     1              0            0     0
========================================================================================================================
Average time to get PetscTime(): 1.19209e-07
Average time for MPI_Barrier(): 4.81606e-06
Average time for zero size MPI_Send(): 1.35899e-05
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
| Time:           Thu Jan 31 22:16:52 2013                                                                             |
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
 --------------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=8.04837, Active time=7.9498                                                        |
 --------------------------------------------------------------------------------------------------------------------
| Event                                  nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                                  w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|--------------------------------------------------------------------------------------------------------------------|
|                                                                                                                    |
|                                                                                                                    |
| CondensedEigenSystem                                                                                               |
|   get_eigenpair()                      14        0.0084      0.000600    0.5441      0.038861    0.11     6.84     |
|   solve()                              14        0.0287      0.002048    5.5454      0.396099    0.36     69.76    |
|                                                                                                                    |
| DofMap                                                                                                             |
|   add_neighbors_to_send_list()         2         0.0062      0.003087    0.0100      0.004992    0.08     0.13     |
|   build_constraint_matrix()            3060      0.0193      0.000006    0.0193      0.000006    0.24     0.24     |
|   build_sparsity()                     2         0.0040      0.001994    0.0124      0.006200    0.05     0.16     |
|   cnstrn_elem_mat_vec()                3060      0.0228      0.000007    0.0228      0.000007    0.29     0.29     |
|   create_dof_constraints()             2         0.0121      0.006073    0.0466      0.023315    0.15     0.59     |
|   distribute_dofs()                    2         0.0222      0.011094    0.0724      0.036203    0.28     0.91     |
|   dof_indices()                        6917      0.3661      0.000053    0.3661      0.000053    4.61     4.61     |
|   prepare_send_list()                  2         0.0001      0.000033    0.0001      0.000033    0.00     0.00     |
|   reinit()                             2         0.0473      0.023629    0.0473      0.023629    0.59     0.59     |
|                                                                                                                    |
| FE                                                                                                                 |
|   compute_shape_functions()            7920      0.1944      0.000025    0.1944      0.000025    2.45     2.45     |
|   init_shape_functions()               1920      0.0603      0.000031    0.0603      0.000031    0.76     0.76     |
|   inverse_map()                        3600      0.0342      0.000010    0.0342      0.000010    0.43     0.43     |
|                                                                                                                    |
| FEMap                                                                                                              |
|   compute_affine_map()                 7920      0.0780      0.000010    0.0780      0.000010    0.98     0.98     |
|   compute_face_map()                   1800      0.0344      0.000019    0.0695      0.000039    0.43     0.87     |
|   init_face_shape_functions()          120       0.0019      0.000016    0.0019      0.000016    0.02     0.02     |
|   init_reference_to_physical_map()     1920      0.0315      0.000016    0.0315      0.000016    0.40     0.40     |
|                                                                                                                    |
| Mesh                                                                                                               |
|   find_neighbors()                     1         0.0148      0.014804    0.0149      0.014934    0.19     0.19     |
|   renumber_nodes_and_elem()            2         0.0006      0.000286    0.0006      0.000286    0.01     0.01     |
|                                                                                                                    |
| MeshCommunication                                                                                                  |
|   assign_global_indices()              1         0.0274      0.027429    0.0297      0.029732    0.35     0.37     |
|   compute_hilbert_indices()            2         0.0108      0.005409    0.0108      0.005409    0.14     0.14     |
|   find_global_indices()                2         0.0043      0.002164    0.0201      0.010034    0.05     0.25     |
|   parallel_sort()                      2         0.0031      0.001525    0.0036      0.001812    0.04     0.05     |
|                                                                                                                    |
| MeshTools::Generation                                                                                              |
|   build_cube()                         1         0.0032      0.003248    0.0032      0.003248    0.04     0.04     |
|                                                                                                                    |
| MetisPartitioner                                                                                                   |
|   partition()                          1         0.0495      0.049516    0.0587      0.058728    0.62     0.74     |
|                                                                                                                    |
| Parallel                                                                                                           |
|   allgather()                          17        0.0008      0.000050    0.0009      0.000054    0.01     0.01     |
|   barrier()                            1         0.0000      0.000021    0.0000      0.000021    0.00     0.00     |
|   broadcast()                          36        0.0004      0.000012    0.0004      0.000012    0.01     0.01     |
|   max(bool)                            7         0.0000      0.000006    0.0000      0.000006    0.00     0.00     |
|   max(scalar)                          168       0.0014      0.000008    0.0014      0.000008    0.02     0.02     |
|   max(vector)                          37        0.0005      0.000014    0.0013      0.000035    0.01     0.02     |
|   maxloc(scalar)                       17        0.0025      0.000150    0.0025      0.000150    0.03     0.03     |
|   min(bool)                            183       0.0012      0.000007    0.0012      0.000007    0.02     0.02     |
|   min(scalar)                          149       0.0102      0.000068    0.0102      0.000068    0.13     0.13     |
|   min(vector)                          37        0.0006      0.000016    0.0014      0.000039    0.01     0.02     |
|   probe()                              338       0.0020      0.000006    0.0020      0.000006    0.02     0.02     |
|   receive()                            290       0.0019      0.000006    0.0037      0.000013    0.02     0.05     |
|   send()                               246       0.0008      0.000003    0.0008      0.000003    0.01     0.01     |
|   send_receive()                       250       0.0019      0.000008    0.0066      0.000026    0.02     0.08     |
|   sum()                                57        0.5364      0.009411    0.5373      0.009427    6.75     6.76     |
|                                                                                                                    |
| Parallel::Request                                                                                                  |
|   wait()                               246       0.0005      0.000002    0.0005      0.000002    0.01     0.01     |
|                                                                                                                    |
| Partitioner                                                                                                        |
|   set_node_processor_ids()             1         0.0016      0.001582    0.0022      0.002206    0.02     0.03     |
|   set_parent_processor_ids()           1         0.0011      0.001147    0.0011      0.001147    0.01     0.01     |
|                                                                                                                    |
| PetscLinearSolver                                                                                                  |
|   solve()                              37        0.1868      0.005048    0.1868      0.005048    2.35     2.35     |
|                                                                                                                    |
| RBConstruction                                                                                                     |
|   add_scaled_Aq()                      51        0.0055      0.000107    1.0737      0.021052    0.07     13.51    |
|   add_scaled_matrix_and_vector()       60        0.4358      0.007263    1.2446      0.020743    5.48     15.66    |
|   clear()                              1         0.0010      0.001010    0.0010      0.001010    0.01     0.01     |
|   compute_Fq_representor_innerprods()  1         0.0008      0.000788    0.0048      0.004806    0.01     0.06     |
|   compute_max_error_bound()            9         0.0025      0.000283    0.0342      0.003801    0.03     0.43     |
|   compute_output_dual_innerprods()     1         0.0021      0.002134    0.0187      0.018676    0.03     0.23     |
|   enrich_RB_space()                    8         0.0046      0.000572    0.0046      0.000572    0.06     0.06     |
|   train_reduced_basis()                1         0.0034      0.003389    0.2963      0.296317    0.04     3.73     |
|   truth_assembly()                     8         0.0168      0.002100    0.0168      0.002100    0.21     0.21     |
|   truth_solve()                        8         0.0015      0.000187    0.0612      0.007645    0.02     0.77     |
|   update_RB_system_matrices()          8         0.0121      0.001508    0.0121      0.001508    0.15     0.15     |
|   update_residual_terms()              8         0.0337      0.004208    0.1570      0.019629    0.42     1.98     |
|                                                                                                                    |
| RBEvaluation                                                                                                       |
|   clear()                              1         0.0001      0.000124    0.0001      0.000124    0.00     0.00     |
|   compute_residual_dual_norm()         81        0.0156      0.000192    0.0156      0.000192    0.20     0.20     |
|   rb_solve()                           81        0.0048      0.000059    0.0303      0.000374    0.06     0.38     |
|   resize_data_structures()             1         0.0004      0.000386    0.0004      0.000386    0.00     0.00     |
|   write_offline_data_to_files()        1         0.0013      0.001284    0.0013      0.001284    0.02     0.02     |
|   write_out_basis_functions()          1         0.0000      0.000050    0.0496      0.049580    0.00     0.62     |
|   write_out_vectors()                  1         0.0188      0.018815    0.0495      0.049530    0.24     0.62     |
|                                                                                                                    |
| RBSCMConstruction                                                                                                  |
|   add_scaled_symm_Aq()                 51        0.0003      0.000007    1.0740      0.021059    0.00     13.51    |
|   compute_SCM_bounding_box()           1         0.0008      0.000845    2.3490      2.349034    0.01     29.55    |
|   compute_SCM_bounds_on_training_set() 8         0.0021      0.000261    0.0117      0.001464    0.03     0.15     |
|   enrich_C_J()                         8         0.0011      0.000135    0.0015      0.000183    0.01     0.02     |
|   evaluate_stability_constant()        8         0.0113      0.001415    4.8266      0.603325    0.14     60.71    |
|   perform_SCM_greedy()                 1         0.0024      0.002404    7.1912      7.191219    0.03     90.46    |
|                                                                                                                    |
| RBSCMEvaluation                                                                                                    |
|   get_SCM_LB()                         153       0.0165      0.000108    0.0165      0.000108    0.21     0.21     |
|   get_SCM_UB()                         72        0.0013      0.000018    0.0013      0.000018    0.02     0.02     |
|   write_offline_data_to_files()        1         0.0002      0.000244    0.0002      0.000244    0.00     0.00     |
|                                                                                                                    |
| SlepcEigenSolver                                                                                                   |
|   solve_generalized()                  14        5.5167      0.394050    5.5167      0.394050    69.39    69.39    |
 --------------------------------------------------------------------------------------------------------------------
| Totals:                                41055     7.9498                                          100.00            |
 --------------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example reduced_basis_ex2:
*  mpirun -np 12 example-devel -online_mode 0 -eps_type lapack -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
***************************************************************
***************************************************************
* Running Example reduced_basis_ex2:
*  mpirun -np 12 example-devel -online_mode 1 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
***************************************************************
 
*** Warning, This code is untested, experimental, or likely to see future API changes: src/reduced_basis/rb_parametrized.C, line 41, compiled Jan 31 2013 at 21:51:32 ***
 EquationSystems
  n_systems()=2
   System #0, "RBConvectionDiffusion"
    Type "RBConstruction"
    Variables="u" 
    Finite Element Types="LAGRANGE", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" 
    Approximation Orders="FIRST", "THIRD" 
    n_dofs()=676
    n_local_dofs()=68
    n_constrained_dofs()=26
    n_local_constrained_dofs()=6
    n_vectors()=1
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 7.99852
      Average Off-Processor Bandwidth <= 1.20118
      Maximum  On-Processor Bandwidth <= 11
      Maximum Off-Processor Bandwidth <= 9
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
    n_local_dofs()=68
    n_constrained_dofs()=0
    n_local_constrained_dofs()=0
    n_vectors()=0
    n_matrices()=2
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 7.99852
      Average Off-Processor Bandwidth <= 1.20118
      Maximum  On-Processor Bandwidth <= 11
      Maximum Off-Processor Bandwidth <= 9
    DofMap Constraints
      Number of DoF Constraints = 0
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

mu_0: 2.000000e-01
mu_1: 7.000000e-01
mu_2: 1.000000e-01

output 1, value = 2.35241, bound = 0.00116972
output 2, value = 0.944951, bound = 0.000720172
output 3, value = 0.944951, bound = 0.000720172
output 4, value = 2.35241, bound = 0.00116972

************************************************************************************************************************
***             WIDEN YOUR WINDOW TO 120 CHARACTERS.  Use 'enscript -r -fCourier9' to print this document            ***
************************************************************************************************************************

---------------------------------------------- PETSc Performance Summary: ----------------------------------------------

/workspace/libmesh/examples/reduced_basis/reduced_basis_ex2/.libs/lt-example-devel on a intel-12. named hbar.ices.utexas.edu with 12 processors, by benkirk Thu Jan 31 22:16:53 2013
Using Petsc Release Version 3.3.0, Patch 2, Fri Jul 13 15:42:00 CDT 2012 

                         Max       Max/Min        Avg      Total 
Time (sec):           3.515e-01      1.00001   3.515e-01
Objects:              3.500e+01      1.00000   3.500e+01
Flops:                8.160e+02      1.78947   6.760e+02  8.112e+03
Flops/sec:            2.321e+03      1.78947   1.923e+03  2.308e+04
MPI Messages:         3.900e+01      3.25000   2.500e+01  3.000e+02
MPI Message Lengths:  1.548e+03      2.37423   4.165e+01  1.250e+04
MPI Reductions:       8.100e+01      1.00000

Flop counting convention: 1 flop = 1 real number operation of type (multiply/divide/add/subtract)
                            e.g., VecAXPY() for real vectors of length N --> 2N flops
                            and VecAXPY() for complex vectors of length N --> 8N flops

Summary of Stages:   ----- Time ------  ----- Flops -----  --- Messages ---  -- Message Lengths --  -- Reductions --
                        Avg     %Total     Avg     %Total   counts   %Total     Avg         %Total   counts   %Total 
 0:      Main Stage: 3.5147e-01 100.0%  8.1120e+03 100.0%  3.000e+02 100.0%  4.165e+01      100.0%  8.000e+01  98.8% 

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

VecCopy                3 1.0 2.3127e-05 2.8 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecSet                16 1.0 1.5736e-05 2.6 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecAXPY                6 1.0 4.9829e-05 1.9 8.16e+02 1.8 0.0e+00 0.0e+00 0.0e+00  0100  0  0  0   0100  0  0  0   163
VecAssemblyBegin      11 1.0 7.1614e-0323.8 0.00e+00 0.0 0.0e+00 0.0e+00 3.3e+01  2  0  0  0 41   2  0  0  0 41     0
VecAssemblyEnd        11 1.0 3.8147e-05 1.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecScatterBegin        2 1.0 8.4877e-05 1.2 0.00e+00 0.0 1.0e+02 8.1e+01 0.0e+00  0  0 33 65  0   0  0 33 65  0     0
VecScatterEnd          2 1.0 2.8133e-05 1.7 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatZeroEntries         6 1.0 2.0504e-05 2.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
------------------------------------------------------------------------------------------------------------------------

Memory usage is given in bytes:

Object Type          Creations   Destructions     Memory  Descendants' Mem.
Reports information only for process 0.

--- Event Stage 0: Main Stage

              Vector    17             17        33784     0
      Vector Scatter     2              2         2072     0
           Index Set     4              4         3120     0
   IS L to G Mapping     2              2         1128     0
              Matrix     9              9        49776     0
              Viewer     1              0            0     0
========================================================================================================================
Average time to get PetscTime(): 1.19209e-07
Average time for MPI_Barrier(): 3.19481e-06
Average time for zero size MPI_Send(): 1.35899e-05
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
| Time:           Thu Jan 31 22:16:53 2013                                                                             |
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
| libMesh Performance: Alive time=0.478045, Active time=0.325008                                               |
 --------------------------------------------------------------------------------------------------------------
| Event                            nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                            w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|--------------------------------------------------------------------------------------------------------------|
|                                                                                                              |
|                                                                                                              |
| DofMap                                                                                                       |
|   add_neighbors_to_send_list()   2         0.0061      0.003060    0.0099      0.004935    1.88     3.04     |
|   build_sparsity()               2         0.0037      0.001871    0.0122      0.006115    1.15     3.76     |
|   create_dof_constraints()       2         0.0118      0.005895    0.0460      0.023005    3.63     14.16    |
|   distribute_dofs()              2         0.0219      0.010944    0.0727      0.036334    6.73     22.36    |
|   dof_indices()                  1001      0.0541      0.000054    0.0541      0.000054    16.65    16.65    |
|   prepare_send_list()            2         0.0001      0.000033    0.0001      0.000033    0.02     0.02     |
|   reinit()                       2         0.0466      0.023293    0.0466      0.023293    14.33    14.33    |
|                                                                                                              |
| EquationSystems                                                                                              |
|   build_solution_vector()        2         0.0027      0.001372    0.0155      0.007745    0.84     4.77     |
|                                                                                                              |
| ExodusII_IO                                                                                                  |
|   write_nodal_data()             2         0.0128      0.006424    0.0128      0.006424    3.95     3.95     |
|                                                                                                              |
| Mesh                                                                                                         |
|   find_neighbors()               1         0.0136      0.013598    0.0142      0.014238    4.18     4.38     |
|   renumber_nodes_and_elem()      2         0.0006      0.000292    0.0006      0.000292    0.18     0.18     |
|                                                                                                              |
| MeshCommunication                                                                                            |
|   assign_global_indices()        1         0.0273      0.027285    0.0295      0.029497    8.40     9.08     |
|   compute_hilbert_indices()      2         0.0120      0.005988    0.0120      0.005988    3.69     3.69     |
|   find_global_indices()          2         0.0043      0.002130    0.0213      0.010659    1.31     6.56     |
|   parallel_sort()                2         0.0031      0.001568    0.0037      0.001826    0.96     1.12     |
|                                                                                                              |
| MeshOutput                                                                                                   |
|   write_equation_systems()       2         0.0002      0.000079    0.0287      0.014365    0.05     8.84     |
|                                                                                                              |
| MeshTools::Generation                                                                                        |
|   build_cube()                   1         0.0031      0.003108    0.0031      0.003108    0.96     0.96     |
|                                                                                                              |
| MetisPartitioner                                                                                             |
|   partition()                    1         0.0490      0.049038    0.0583      0.058339    15.09    17.95    |
|                                                                                                              |
| Parallel                                                                                                     |
|   allgather()                    17        0.0019      0.000110    0.0019      0.000114    0.58     0.60     |
|   barrier()                      1         0.0000      0.000031    0.0000      0.000031    0.01     0.01     |
|   broadcast()                    13        0.0003      0.000025    0.0003      0.000021    0.10     0.08     |
|   max(bool)                      3         0.0000      0.000006    0.0000      0.000006    0.01     0.01     |
|   max(scalar)                    186       0.0012      0.000006    0.0012      0.000006    0.37     0.37     |
|   max(vector)                    41        0.0005      0.000013    0.0013      0.000031    0.16     0.40     |
|   min(bool)                      213       0.0013      0.000006    0.0013      0.000006    0.40     0.40     |
|   min(scalar)                    175       0.0100      0.000057    0.0100      0.000057    3.07     3.07     |
|   min(vector)                    41        0.0006      0.000016    0.0015      0.000037    0.20     0.46     |
|   probe()                        290       0.0023      0.000008    0.0023      0.000008    0.71     0.71     |
|   receive()                      268       0.0017      0.000006    0.0040      0.000015    0.52     1.22     |
|   send()                         268       0.0008      0.000003    0.0008      0.000003    0.26     0.26     |
|   send_receive()                 250       0.0020      0.000008    0.0070      0.000028    0.60     2.16     |
|   sum()                          44        0.0009      0.000020    0.0023      0.000052    0.27     0.71     |
|                                                                                                              |
| Parallel::Request                                                                                            |
|   wait()                         268       0.0005      0.000002    0.0005      0.000002    0.16     0.16     |
|                                                                                                              |
| Partitioner                                                                                                  |
|   set_node_processor_ids()       1         0.0015      0.001538    0.0028      0.002828    0.47     0.87     |
|   set_parent_processor_ids()     1         0.0011      0.001147    0.0011      0.001147    0.35     0.35     |
|                                                                                                              |
| RBConstruction                                                                                               |
|   clear()                        1         0.0004      0.000405    0.0004      0.000405    0.12     0.12     |
|   load_basis_function()          1         0.0001      0.000092    0.0001      0.000092    0.03     0.03     |
|   load_rb_solution()             1         0.0003      0.000350    0.0003      0.000350    0.11     0.11     |
|                                                                                                              |
| RBEvaluation                                                                                                 |
|   clear()                        1         0.0001      0.000118    0.0001      0.000118    0.04     0.04     |
|   compute_residual_dual_norm()   1         0.0003      0.000311    0.0003      0.000311    0.10     0.10     |
|   rb_solve()                     1         0.0106      0.010619    0.0114      0.011374    3.27     3.50     |
|   read_in_basis_functions()      1         0.0000      0.000046    0.0414      0.041363    0.01     12.73    |
|   read_in_vectors()              1         0.0108      0.010799    0.0413      0.041316    3.32     12.71    |
|   read_offline_data_from_files() 1         0.0012      0.001197    0.0018      0.001788    0.37     0.55     |
|   resize_data_structures()       1         0.0006      0.000590    0.0006      0.000590    0.18     0.18     |
|                                                                                                              |
| RBSCMEvaluation                                                                                              |
|   get_SCM_LB()                   1         0.0004      0.000444    0.0004      0.000444    0.14     0.14     |
|   read_offline_data_from_files() 1         0.0003      0.000270    0.0003      0.000270    0.08     0.08     |
 --------------------------------------------------------------------------------------------------------------
| Totals:                          3123      0.3250                                          100.00            |
 --------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example reduced_basis_ex2:
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
