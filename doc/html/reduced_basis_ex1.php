<?php $root=""; ?>
<?php require($root."navigation.php"); ?>
<html>
<head>
  <?php load_style($root); ?>
</head>
 
<body>
 
<?php make_navigation("reduced_basis_ex1",$root)?>
 
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
        struct CDRBThetaExpansion : RBThetaExpansion
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
        struct CDRBAssemblyExpansion : RBAssemblyExpansion
        {
        
          /**
           * Constructor.
           */
          CDRBAssemblyExpansion()
            :
            L0(0.7,0.8,0.7,0.8),
            L1(0.2,0.3,0.7,0.8),
            L2(0.2,0.3,0.2,0.3),
            L3(0.7,0.8,0.2,0.3)
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
        class SimpleRBConstruction : public RBConstruction
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
          virtual ~SimpleRBConstruction () { }
        
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
The theta expansion comes from the RBEvaluation object.
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
<br><br><br> <h1> The source file reduced_basis_ex1.C with comments: </h1> 
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
<h1>Reduced Basis Example 1 - Certified Reduced Basis Method</h1>
 

<br><br>In this example problem we use the Certified Reduced Basis method
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
specified in reduced_basis_ex1.in.


<br><br>We also attach four outputs to the system which are averages over certain
subregions of the domain. In Online mode, we print out the values of these
outputs as well as rigorous error bounds with respect to the output
associated with the "truth" finite element discretization.


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
Parse the input file (reduced_basis_ex1.in) using GetPot
</div>

<div class ="fragment">
<pre>
          std::string parameters_filename = "reduced_basis_ex1.in";
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
  
  <B><FONT COLOR="#228B22">struct</FONT></B> ThetaA0 : RBTheta { <B><FONT COLOR="#228B22">virtual</FONT></B> Number evaluate(<B><FONT COLOR="#228B22">const</FONT></B> RBParameters&amp; )   { <B><FONT COLOR="#A020F0">return</FONT></B> 0.05;  } };
  <B><FONT COLOR="#228B22">struct</FONT></B> ThetaA1 : RBTheta { <B><FONT COLOR="#228B22">virtual</FONT></B> Number evaluate(<B><FONT COLOR="#228B22">const</FONT></B> RBParameters&amp; mu) { <B><FONT COLOR="#A020F0">return</FONT></B> mu.get_value(<B><FONT COLOR="#BC8F8F">&quot;x_vel&quot;</FONT></B>); } };
  <B><FONT COLOR="#228B22">struct</FONT></B> ThetaA2 : RBTheta { <B><FONT COLOR="#228B22">virtual</FONT></B> Number evaluate(<B><FONT COLOR="#228B22">const</FONT></B> RBParameters&amp; mu) { <B><FONT COLOR="#A020F0">return</FONT></B> mu.get_value(<B><FONT COLOR="#BC8F8F">&quot;y_vel&quot;</FONT></B>); } };
  
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
  
  <B><FONT COLOR="#228B22">struct</FONT></B> CDRBThetaExpansion : RBThetaExpansion
  {
  
    <I><FONT COLOR="#B22222">/**
     * Constructor.
     */</FONT></I>
    CDRBThetaExpansion()
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
  
  <B><FONT COLOR="#228B22">struct</FONT></B> CDRBAssemblyExpansion : RBAssemblyExpansion
  {
  
    <I><FONT COLOR="#B22222">/**
     * Constructor.
     */</FONT></I>
    CDRBAssemblyExpansion()
      :
      L0(0.7,0.8,0.7,0.8),
      L1(0.2,0.3,0.7,0.8),
      L2(0.2,0.3,0.2,0.3),
      L3(0.7,0.8,0.2,0.3)
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
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/rb_construction.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/fe_base.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;assembly.h&quot;</FONT></B>
  
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
  
  <B><FONT COLOR="#228B22">class</FONT></B> SimpleRBConstruction : <B><FONT COLOR="#228B22">public</FONT></B> RBConstruction
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
    <B><FONT COLOR="#228B22">virtual</FONT></B> ~SimpleRBConstruction () { }
  
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
      
      dirichlet_bc-&gt;b.insert(0);
      dirichlet_bc-&gt;b.insert(1);
      dirichlet_bc-&gt;b.insert(2);
      dirichlet_bc-&gt;b.insert(3);
      dirichlet_bc-&gt;variables.push_back(u_var);
      
      get_dof_map().add_dirichlet_boundary(*dirichlet_bc);
  
      <B><FONT COLOR="#5F9EA0">Parent</FONT></B>::init_data();
  
      set_rb_assembly_expansion(cd_rb_assembly_expansion);
  
      set_inner_product_assembly(cd_rb_assembly_expansion.A0_assembly);
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
<br><br><br> <h1> The source file reduced_basis_ex1.C without comments: </h1> 
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
  
  #<B><FONT COLOR="#A020F0">if</FONT></B> !defined(LIBMESH_HAVE_XDR)
    libmesh_example_assert(false, <B><FONT COLOR="#BC8F8F">&quot;--enable-xdr&quot;</FONT></B>);
  #elif defined(LIBMESH_DEFAULT_SINGLE_PRECISION)
    libmesh_example_assert(false, <B><FONT COLOR="#BC8F8F">&quot;--disable-singleprecision&quot;</FONT></B>);
  #endif
    libmesh_example_assert(libMesh::default_solver_package() == PETSC_SOLVERS, <B><FONT COLOR="#BC8F8F">&quot;--enable-petsc&quot;</FONT></B>);
  
    libmesh_example_assert(2 &lt;= LIBMESH_DIM, <B><FONT COLOR="#BC8F8F">&quot;2D support&quot;</FONT></B>);
    
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::string parameters_filename = <B><FONT COLOR="#BC8F8F">&quot;reduced_basis_ex1.in&quot;</FONT></B>;
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
  }
  
</pre> 
<a name="output"></a> 
<br><br><br> <h1> The console output of the program: </h1> 
<pre>
***************************************************************
* Running Example reduced_basis_ex1:
*  mpirun -np 12 example-devel -online_mode 0 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
***************************************************************
 
*** Warning, This code is untested, experimental, or likely to see future API changes: src/reduced_basis/rb_parametrized.C, line 41, compiled Jan 31 2013 at 21:51:32 ***
 EquationSystems
  n_systems()=1
   System #0, "RBConvectionDiffusion"
    Type "RBConstruction"
    Variables="u" 
    Finite Element Types="LAGRANGE", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" 
    Approximation Orders="FIRST", "THIRD" 
    n_dofs()=2601
    n_local_dofs()=247
    n_constrained_dofs()=200
    n_local_constrained_dofs()=34
    n_vectors()=1
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 8.47174
      Average Off-Processor Bandwidth <= 0.625144
      Maximum  On-Processor Bandwidth <= 11
      Maximum Off-Processor Bandwidth <= 6
    DofMap Constraints
      Number of DoF Constraints = 200
      Average DoF Constraint Length= 0
      Number of Node Constraints = 0

 Mesh Information:
  mesh_dimension()=2
  spatial_dimension()=3
  n_nodes()=2601
    n_local_nodes()=247
  n_elem()=2500
    n_local_elem()=210
    n_active_elem()=2500
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
Parameter x_vel: Min = -2, Max = 2, value = -2
Parameter y_vel: Min = -2, Max = 2, value = 0.79
n_training_samples: 100
single-matrix mode? 0
reuse preconditioner? 1
use a relative error bound in greedy? 0
write out data during basis training? 0
quiet mode? 1

Compute output dual inner products
output_dual_innerprods[0][0] = 0.323675
output_dual_innerprods[1][0] = 0.323675
output_dual_innerprods[2][0] = 0.323675
output_dual_innerprods[3][0] = 0.323675

---- Performing Greedy basis enrichment ----

---- Basis dimension: 0 ----
Performing RB solves on training set
Maximum (absolute) error bound is 3.74824

Performing truth solve at parameter:
x_vel: -2.000000e+00
y_vel: -2.000000e+00

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 1 ----
Performing RB solves on training set
Maximum (absolute) error bound is 6.67625

Performing truth solve at parameter:
x_vel: 2.000000e+00
y_vel: 2.000000e+00

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 2 ----
Performing RB solves on training set
Maximum (absolute) error bound is 5.18295

Performing truth solve at parameter:
x_vel: 2.000000e+00
y_vel: -2.000000e+00

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 3 ----
Performing RB solves on training set
Maximum (absolute) error bound is 3.59028

Performing truth solve at parameter:
x_vel: -2.000000e+00
y_vel: 2.000000e+00

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 4 ----
Performing RB solves on training set
Maximum (absolute) error bound is 2.71338

Performing truth solve at parameter:
x_vel: 2.222222e-01
y_vel: 2.222222e-01

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 5 ----
Performing RB solves on training set
Maximum (absolute) error bound is 2.50784

Performing truth solve at parameter:
x_vel: -2.222222e-01
y_vel: -6.666667e-01

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 6 ----
Performing RB solves on training set
Maximum (absolute) error bound is 2.15653

Performing truth solve at parameter:
x_vel: -6.666667e-01
y_vel: 2.222222e-01

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 7 ----
Performing RB solves on training set
Maximum (absolute) error bound is 1.71352

Performing truth solve at parameter:
x_vel: 1.111111e+00
y_vel: -2.222222e-01

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 8 ----
Performing RB solves on training set
Maximum (absolute) error bound is 1.61442

Performing truth solve at parameter:
x_vel: -2.222222e-01
y_vel: 1.555556e+00

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 9 ----
Performing RB solves on training set
Maximum (absolute) error bound is 1.09712

Performing truth solve at parameter:
x_vel: 2.222222e-01
y_vel: -2.000000e+00

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 10 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.9794

Performing truth solve at parameter:
x_vel: -2.000000e+00
y_vel: -2.222222e-01

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 11 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.763958

Performing truth solve at parameter:
x_vel: 2.222222e-01
y_vel: -2.222222e-01

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 12 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.616746

Performing truth solve at parameter:
x_vel: 1.111111e+00
y_vel: 1.111111e+00

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 13 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.530181

Performing truth solve at parameter:
x_vel: 2.000000e+00
y_vel: 2.222222e-01

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 14 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.475932

Performing truth solve at parameter:
x_vel: -2.222222e-01
y_vel: 6.666667e-01

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 15 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.440399

Performing truth solve at parameter:
x_vel: 6.666667e-01
y_vel: -1.111111e+00

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 16 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.423537

Performing truth solve at parameter:
x_vel: -1.111111e+00
y_vel: 1.111111e+00

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 17 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.408558

Performing truth solve at parameter:
x_vel: -1.111111e+00
y_vel: -6.666667e-01

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 18 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.333764

Performing truth solve at parameter:
x_vel: -2.222222e-01
y_vel: -2.222222e-01

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 19 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.196867

Performing truth solve at parameter:
x_vel: 6.666667e-01
y_vel: 2.222222e-01

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 20 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.191475

Maximum number of basis functions reached: Nmax = 20
************************************************************************************************************************
***             WIDEN YOUR WINDOW TO 120 CHARACTERS.  Use 'enscript -r -fCourier9' to print this document            ***
************************************************************************************************************************

---------------------------------------------- PETSc Performance Summary: ----------------------------------------------

/workspace/libmesh/examples/reduced_basis/reduced_basis_ex1/.libs/lt-example-devel on a intel-12. named hbar.ices.utexas.edu with 12 processors, by benkirk Thu Jan 31 22:16:14 2013
Using Petsc Release Version 3.3.0, Patch 2, Fri Jul 13 15:42:00 CDT 2012 

                         Max       Max/Min        Avg      Total 
Time (sec):           2.597e+00      1.00000   2.597e+00
Objects:              1.148e+03      1.00000   1.148e+03
Flops:                2.136e+08      1.37641   1.926e+08  2.311e+09
Flops/sec:            8.224e+07      1.37641   7.416e+07  8.899e+08
MPI Messages:         7.785e+04      3.50000   4.448e+04  5.338e+05
MPI Message Lengths:  6.110e+06      2.19094   9.951e+01  5.312e+07
MPI Reductions:       2.959e+04      1.00000

Flop counting convention: 1 flop = 1 real number operation of type (multiply/divide/add/subtract)
                            e.g., VecAXPY() for real vectors of length N --> 2N flops
                            and VecAXPY() for complex vectors of length N --> 8N flops

Summary of Stages:   ----- Time ------  ----- Flops -----  --- Messages ---  -- Message Lengths --  -- Reductions --
                        Avg     %Total     Avg     %Total   counts   %Total     Avg         %Total   counts   %Total 
 0:      Main Stage: 2.5971e+00 100.0%  2.3113e+09 100.0%  5.338e+05 100.0%  9.951e+01      100.0%  2.959e+04 100.0% 

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

VecDot              4075 1.0 4.1675e-02 1.6 2.01e+06 1.4 0.0e+00 0.0e+00 4.1e+03  1  1  0  0 14   1  1  0  0 14   507
VecMDot             6694 1.0 1.2057e-01 1.4 4.85e+07 1.4 0.0e+00 0.0e+00 6.7e+03  4 22  0  0 23   4 22  0  0 23  4232
VecNorm             7021 1.0 6.9268e-02 1.1 3.47e+06 1.4 0.0e+00 0.0e+00 7.0e+03  3  2  0  0 24   3  2  0  0 24   527
VecScale            7036 1.0 4.4811e-03 1.5 1.73e+06 1.4 0.0e+00 0.0e+00 0.0e+00  0  1  0  0  0   0  1  0  0  0  4072
VecCopy              622 1.0 5.3024e-04 1.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecSet             12793 1.0 6.1331e-03 1.3 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecAXPY              699 1.0 9.8090e-0317.3 3.45e+05 1.4 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0   371
VecMAXPY            6936 1.0 3.1447e-02 1.4 5.19e+07 1.4 0.0e+00 0.0e+00 0.0e+00  1 24  0  0  0   1 24  0  0  0 17372
VecAssemblyBegin     595 1.0 8.6834e-02 5.3 0.00e+00 0.0 2.4e+02 1.6e+02 1.8e+03  2  0  0  0  6   2  0  0  0  6     0
VecAssemblyEnd       595 1.0 4.1270e-04 1.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecScatterBegin    10936 1.0 4.7227e-02 1.9 0.00e+00 0.0 5.2e+05 1.0e+02 0.0e+00  1  0 98 99  0   1  0 98 99  0     0
VecScatterEnd      10936 1.0 5.3661e-02 1.8 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  2  0  0  0  0   2  0  0  0  0     0
VecNormalize        6936 1.0 7.4718e-02 1.1 5.14e+06 1.4 0.0e+00 0.0e+00 6.9e+03  3  2  0  0 23   3  2  0  0 23   724
MatMult             6936 1.0 1.0532e-01 1.3 2.77e+07 1.4 3.3e+05 1.0e+02 0.0e+00  4 13 62 63  0   4 13 62 63  0  2832
MatMultAdd          3915 1.0 5.1466e-02 1.4 1.66e+07 1.4 1.9e+05 1.0e+02 0.0e+00  2  8 35 35  0   2  8 35 35  0  3469
MatSolve            7021 1.0 6.1690e-02 1.3 6.26e+07 1.4 0.0e+00 0.0e+00 0.0e+00  2 29  0  0  0   2 29  0  0  0 11015
MatLUFactorNum        42 1.0 5.2474e-03 1.4 1.86e+06 1.7 0.0e+00 0.0e+00 0.0e+00  0  1  0  0  0   0  1  0  0  0  3594
MatILUFactorSym       42 1.0 1.0757e-02 1.4 0.00e+00 0.0 0.0e+00 0.0e+00 1.3e+02  0  0  0  0  0   0  0  0  0  0     0
MatAssemblyBegin    4213 1.0 1.4438e-01 1.8 0.00e+00 0.0 2.9e+02 5.0e+02 8.4e+03  4  0  0  0 28   4  0  0  0 28     0
MatAssemblyEnd      4213 1.0 4.6914e-02 1.2 0.00e+00 0.0 8.3e+03 2.7e+01 7.0e+02  2  0  2  0  2   2  0  2  0  2     0
MatGetRow          40508 1.4 6.7892e-03 1.4 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatGetRowIJ           42 1.0 2.2650e-05 3.3 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatGetOrdering        42 1.0 1.1733e-03 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 8.4e+01  0  0  0  0  0   0  0  0  0  0     0
MatZeroEntries        56 1.0 1.7405e-04 1.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatAXPY               82 1.0 7.9064e-02 1.0 0.00e+00 0.0 7.9e+03 2.7e+01 1.3e+03  3  0  1  0  4   3  0  1  0  4     0
KSPGMRESOrthog      6694 1.0 1.5302e-01 1.3 9.70e+07 1.4 0.0e+00 0.0e+00 6.7e+03  5 44  0  0 23   5 44  0  0 23  6677
KSPSetUp             127 1.0 6.7234e-04 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPSolve              85 1.0 4.6043e-01 1.0 1.95e+08 1.4 3.3e+05 1.0e+02 1.4e+04 18 91 62 63 47  18 91 62 63 47  4583
PCSetUp               84 1.0 2.5854e-02 1.2 1.86e+06 1.7 0.0e+00 0.0e+00 2.1e+02  1  1  0  0  1   1  1  0  0  1   729
PCSetUpOnBlocks       85 1.0 2.3565e-02 1.2 1.86e+06 1.7 0.0e+00 0.0e+00 2.1e+02  1  1  0  0  1   1  1  0  0  1   800
PCApply             7021 1.0 1.2615e-01 1.1 6.26e+07 1.4 0.0e+00 0.0e+00 0.0e+00  5 29  0  0  0   5 29  0  0  0  5386
------------------------------------------------------------------------------------------------------------------------

Memory usage is given in bytes:

Object Type          Creations   Destructions     Memory  Descendants' Mem.
Reports information only for process 0.

--- Event Stage 0: Main Stage

              Vector   449            449      1228840     0
      Vector Scatter    88             88        91168     0
           Index Set   302            302       280096     0
   IS L to G Mapping     1              1          564     0
              Matrix   303            303      5903232     0
       Krylov Solver     2              2        19360     0
      Preconditioner     2              2         1784     0
              Viewer     1              0            0     0
========================================================================================================================
Average time to get PetscTime(): 1.19209e-07
Average time for MPI_Barrier(): 4.19617e-06
Average time for zero size MPI_Send(): 1.38283e-05
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
| Time:           Thu Jan 31 22:16:14 2013                                                                             |
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
| libMesh Performance: Alive time=2.66468, Active time=2.53571                                                      |
 -------------------------------------------------------------------------------------------------------------------
| Event                                 nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                                 w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|-------------------------------------------------------------------------------------------------------------------|
|                                                                                                                   |
|                                                                                                                   |
| DofMap                                                                                                            |
|   add_neighbors_to_send_list()        1         0.0122      0.012178    0.0165      0.016472    0.48     0.65     |
|   build_constraint_matrix()           1890      0.0163      0.000009    0.0163      0.000009    0.64     0.64     |
|   build_sparsity()                    1         0.0071      0.007148    0.0208      0.020789    0.28     0.82     |
|   cnstrn_elem_mat_vec()               1890      0.0216      0.000011    0.0216      0.000011    0.85     0.85     |
|   create_dof_constraints()            1         0.0426      0.042559    0.1802      0.180193    1.68     7.11     |
|   distribute_dofs()                   1         0.0416      0.041630    0.1382      0.138202    1.64     5.45     |
|   dof_indices()                       6569      0.3702      0.000056    0.3702      0.000056    14.60    14.60    |
|   prepare_send_list()                 1         0.0001      0.000076    0.0001      0.000076    0.00     0.00     |
|   reinit()                            1         0.0925      0.092471    0.0925      0.092471    3.65     3.65     |
|                                                                                                                   |
| FE                                                                                                                |
|   compute_shape_functions()           4374      0.1165      0.000027    0.1165      0.000027    4.60     4.60     |
|   init_shape_functions()              612       0.0175      0.000029    0.0175      0.000029    0.69     0.69     |
|   inverse_map()                       1188      0.0123      0.000010    0.0123      0.000010    0.48     0.48     |
|                                                                                                                   |
| FEMap                                                                                                             |
|   compute_affine_map()                4374      0.0485      0.000011    0.0485      0.000011    1.91     1.91     |
|   compute_face_map()                  594       0.0118      0.000020    0.0243      0.000041    0.46     0.96     |
|   init_face_shape_functions()         18        0.0003      0.000017    0.0003      0.000017    0.01     0.01     |
|   init_reference_to_physical_map()    612       0.0095      0.000015    0.0095      0.000015    0.37     0.37     |
|                                                                                                                   |
| Mesh                                                                                                              |
|   find_neighbors()                    1         0.0549      0.054886    0.0559      0.055884    2.16     2.20     |
|   renumber_nodes_and_elem()           2         0.0025      0.001275    0.0025      0.001275    0.10     0.10     |
|                                                                                                                   |
| MeshCommunication                                                                                                 |
|   assign_global_indices()             1         0.1027      0.102660    0.1061      0.106099    4.05     4.18     |
|   compute_hilbert_indices()           2         0.0426      0.021302    0.0426      0.021302    1.68     1.68     |
|   find_global_indices()               2         0.0167      0.008339    0.0662      0.033118    0.66     2.61     |
|   parallel_sort()                     2         0.0046      0.002312    0.0056      0.002784    0.18     0.22     |
|                                                                                                                   |
| MeshTools::Generation                                                                                             |
|   build_cube()                        1         0.0107      0.010708    0.0107      0.010708    0.42     0.42     |
|                                                                                                                   |
| MetisPartitioner                                                                                                  |
|   partition()                         1         0.1339      0.133895    0.1660      0.165962    5.28     6.54     |
|                                                                                                                   |
| Parallel                                                                                                          |
|   allgather()                         14        0.0038      0.000271    0.0039      0.000276    0.15     0.15     |
|   barrier()                           1         0.0000      0.000020    0.0000      0.000020    0.00     0.00     |
|   broadcast()                         44        0.0007      0.000017    0.0007      0.000017    0.03     0.03     |
|   max(bool)                           5         0.0000      0.000007    0.0000      0.000007    0.00     0.00     |
|   max(scalar)                         119       0.0010      0.000008    0.0010      0.000008    0.04     0.04     |
|   max(vector)                         28        0.0004      0.000014    0.0010      0.000037    0.02     0.04     |
|   maxloc(scalar)                      21        0.0073      0.000350    0.0073      0.000350    0.29     0.29     |
|   min(bool)                           138       0.0011      0.000008    0.0011      0.000008    0.04     0.04     |
|   min(scalar)                         113       0.0209      0.000185    0.0209      0.000185    0.82     0.82     |
|   min(vector)                         28        0.0005      0.000018    0.0015      0.000053    0.02     0.06     |
|   probe()                             272       0.0025      0.000009    0.0025      0.000009    0.10     0.10     |
|   receive()                           224       0.0020      0.000009    0.0042      0.000019    0.08     0.17     |
|   send()                              180       0.0007      0.000004    0.0007      0.000004    0.03     0.03     |
|   send_receive()                      184       0.0016      0.000009    0.0061      0.000033    0.06     0.24     |
|   sum()                               44        0.0017      0.000039    0.0021      0.000049    0.07     0.08     |
|                                                                                                                   |
| Parallel::Request                                                                                                 |
|   wait()                              180       0.0004      0.000002    0.0004      0.000002    0.01     0.01     |
|                                                                                                                   |
| Partitioner                                                                                                       |
|   set_node_processor_ids()            1         0.0048      0.004814    0.0058      0.005805    0.19     0.23     |
|   set_parent_processor_ids()          1         0.0046      0.004563    0.0046      0.004563    0.18     0.18     |
|                                                                                                                   |
| PetscLinearSolver                                                                                                 |
|   solve()                             85        0.4950      0.005824    0.4950      0.005824    19.52    19.52    |
|                                                                                                                   |
| RBConstruction                                                                                                    |
|   add_scaled_matrix_and_vector()      9         0.1185      0.013165    0.5951      0.066122    4.67     23.47    |
|   clear()                             1         0.0011      0.001128    0.0011      0.001128    0.04     0.04     |
|   compute_Fq_representor_innerprods() 1         0.0011      0.001073    0.0067      0.006689    0.04     0.26     |
|   compute_max_error_bound()           21        0.0051      0.000242    0.2378      0.011324    0.20     9.38     |
|   compute_output_dual_innerprods()    1         0.0037      0.003730    0.0373      0.037295    0.15     1.47     |
|   enrich_RB_space()                   20        0.0194      0.000971    0.0194      0.000971    0.77     0.77     |
|   train_reduced_basis()               1         0.0078      0.007813    1.0469      1.046914    0.31     41.29    |
|   truth_assembly()                    20        0.0712      0.003558    0.0712      0.003558    2.81     2.81     |
|   truth_solve()                       20        0.0038      0.000190    0.2048      0.010240    0.15     8.08     |
|   update_RB_system_matrices()         20        0.0633      0.003166    0.0633      0.003166    2.50     2.50     |
|   update_residual_terms()             20        0.1431      0.007155    0.4691      0.023457    5.64     18.50    |
|                                                                                                                   |
| RBEvaluation                                                                                                      |
|   clear()                             1         0.0003      0.000266    0.0003      0.000266    0.01     0.01     |
|   compute_residual_dual_norm()        189       0.1518      0.000803    0.1518      0.000803    5.99     5.99     |
|   rb_solve()                          189       0.0719      0.000380    0.2238      0.001184    2.83     8.83     |
|   resize_data_structures()            1         0.0006      0.000569    0.0006      0.000569    0.02     0.02     |
|   write_offline_data_to_files()       1         0.0020      0.002024    0.0020      0.002024    0.08     0.08     |
|   write_out_basis_functions()         1         0.0001      0.000071    0.2407      0.240686    0.00     9.49     |
|   write_out_vectors()                 1         0.1327      0.132688    0.2406      0.240615    5.23     9.49     |
 -------------------------------------------------------------------------------------------------------------------
| Totals:                               24338     2.5357                                          100.00            |
 -------------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example reduced_basis_ex1:
*  mpirun -np 12 example-devel -online_mode 0 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
***************************************************************
***************************************************************
* Running Example reduced_basis_ex1:
*  mpirun -np 12 example-devel -online_mode 1 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
***************************************************************
 
*** Warning, This code is untested, experimental, or likely to see future API changes: src/reduced_basis/rb_parametrized.C, line 41, compiled Jan 31 2013 at 21:51:32 ***
 EquationSystems
  n_systems()=1
   System #0, "RBConvectionDiffusion"
    Type "RBConstruction"
    Variables="u" 
    Finite Element Types="LAGRANGE", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" 
    Approximation Orders="FIRST", "THIRD" 
    n_dofs()=2601
    n_local_dofs()=247
    n_constrained_dofs()=200
    n_local_constrained_dofs()=34
    n_vectors()=1
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 8.47174
      Average Off-Processor Bandwidth <= 0.625144
      Maximum  On-Processor Bandwidth <= 11
      Maximum Off-Processor Bandwidth <= 6
    DofMap Constraints
      Number of DoF Constraints = 200
      Average DoF Constraint Length= 0
      Number of Node Constraints = 0

 Mesh Information:
  mesh_dimension()=2
  spatial_dimension()=3
  n_nodes()=2601
    n_local_nodes()=247
  n_elem()=2500
    n_local_elem()=210
    n_active_elem()=2500
  n_subdomains()=1
  n_partitions()=12
  n_processors()=12
  n_threads()=1
  processor_id()=0

x_vel: -2.000000e+00
y_vel: 7.900000e-01

output 1, value = 0.12351, bound = 0.112915
output 2, value = 0.379732, bound = 0.112915
output 3, value = 0.2501, bound = 0.112915
output 4, value = 0.105334, bound = 0.112915

************************************************************************************************************************
***             WIDEN YOUR WINDOW TO 120 CHARACTERS.  Use 'enscript -r -fCourier9' to print this document            ***
************************************************************************************************************************

---------------------------------------------- PETSc Performance Summary: ----------------------------------------------

/workspace/libmesh/examples/reduced_basis/reduced_basis_ex1/.libs/lt-example-devel on a intel-12. named hbar.ices.utexas.edu with 12 processors, by benkirk Thu Jan 31 22:16:15 2013
Using Petsc Release Version 3.3.0, Patch 2, Fri Jul 13 15:42:00 CDT 2012 

                         Max       Max/Min        Avg      Total 
Time (sec):           9.766e-01      1.00819   9.746e-01
Objects:              3.300e+01      1.00000   3.300e+01
Flops:                9.880e+03      1.42775   8.670e+03  1.040e+05
Flops/sec:            1.012e+04      1.41615   8.895e+03  1.067e+05
MPI Messages:         2.800e+01      3.50000   1.600e+01  1.920e+02
MPI Message Lengths:  2.148e+03      2.19632   9.694e+01  1.861e+04
MPI Reductions:       1.070e+02      1.00000

Flop counting convention: 1 flop = 1 real number operation of type (multiply/divide/add/subtract)
                            e.g., VecAXPY() for real vectors of length N --> 2N flops
                            and VecAXPY() for complex vectors of length N --> 8N flops

Summary of Stages:   ----- Time ------  ----- Flops -----  --- Messages ---  -- Message Lengths --  -- Reductions --
                        Avg     %Total     Avg     %Total   counts   %Total     Avg         %Total   counts   %Total 
 0:      Main Stage: 9.7459e-01 100.0%  1.0404e+05 100.0%  1.920e+02 100.0%  9.694e+01      100.0%  1.060e+02  99.1% 

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

VecCopy                3 1.0 2.0981e-05 1.5 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecSet                25 1.0 1.8358e-05 1.5 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecAXPY               20 1.0 5.9605e-05 1.5 9.88e+03 1.4 0.0e+00 0.0e+00 0.0e+00  0100  0  0  0   0100  0  0  0  1746
VecAssemblyBegin      23 1.0 2.9297e-02 3.1 0.00e+00 0.0 0.0e+00 0.0e+00 6.9e+01  3  0  0  0 64   3  0  0  0 65     0
VecAssemblyEnd        23 1.0 5.5075e-05 2.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecScatterBegin        2 1.0 9.2983e-05 1.1 0.00e+00 0.0 9.6e+01 1.5e+02 0.0e+00  0  0 50 79  0   0  0 50 79  0     0
VecScatterEnd          2 1.0 2.6941e-05 1.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatZeroEntries         2 1.0 1.7166e-05 1.7 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
------------------------------------------------------------------------------------------------------------------------

Memory usage is given in bytes:

Object Type          Creations   Destructions     Memory  Descendants' Mem.
Reports information only for process 0.

--- Event Stage 0: Main Stage

              Vector    25             25        85096     0
      Vector Scatter     1              1         1036     0
           Index Set     2              2         1648     0
   IS L to G Mapping     1              1          564     0
              Matrix     3              3        40388     0
              Viewer     1              0            0     0
========================================================================================================================
Average time to get PetscTime(): 1.19209e-07
Average time for MPI_Barrier(): 3.58105e-05
Average time for zero size MPI_Send(): 5.88298e-05
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
| Time:           Thu Jan 31 22:16:15 2013                                                                             |
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
| libMesh Performance: Alive time=1.05739, Active time=0.936553                                                |
 --------------------------------------------------------------------------------------------------------------
| Event                            nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                            w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|--------------------------------------------------------------------------------------------------------------|
|                                                                                                              |
|                                                                                                              |
| DofMap                                                                                                       |
|   add_neighbors_to_send_list()   1         0.0121      0.012061    0.0163      0.016336    1.29     1.74     |
|   build_sparsity()               1         0.0071      0.007080    0.0205      0.020477    0.76     2.19     |
|   create_dof_constraints()       1         0.0404      0.040436    0.1778      0.177804    4.32     18.98    |
|   distribute_dofs()              1         0.0429      0.042917    0.1373      0.137297    4.58     14.66    |
|   dof_indices()                  3209      0.1747      0.000054    0.1747      0.000054    18.66    18.66    |
|   prepare_send_list()            1         0.0001      0.000075    0.0001      0.000075    0.01     0.01     |
|   reinit()                       1         0.0919      0.091878    0.0919      0.091878    9.81     9.81     |
|                                                                                                              |
| EquationSystems                                                                                              |
|   build_solution_vector()        2         0.0059      0.002967    0.0300      0.014976    0.63     3.20     |
|                                                                                                              |
| ExodusII_IO                                                                                                  |
|   write_nodal_data()             2         0.0402      0.020105    0.0402      0.020105    4.29     4.29     |
|                                                                                                              |
| Mesh                                                                                                         |
|   find_neighbors()               1         0.0532      0.053173    0.0591      0.059085    5.68     6.31     |
|   renumber_nodes_and_elem()      2         0.0019      0.000963    0.0019      0.000963    0.21     0.21     |
|                                                                                                              |
| MeshCommunication                                                                                            |
|   assign_global_indices()        1         0.1035      0.103457    0.1056      0.105624    11.05    11.28    |
|   compute_hilbert_indices()      2         0.0424      0.021218    0.0424      0.021218    4.53     4.53     |
|   find_global_indices()          2         0.0163      0.008150    0.0743      0.037150    1.74     7.93     |
|   parallel_sort()                2         0.0048      0.002421    0.0104      0.005206    0.52     1.11     |
|                                                                                                              |
| MeshOutput                                                                                                   |
|   write_equation_systems()       2         0.0002      0.000092    0.0706      0.035298    0.02     7.54     |
|                                                                                                              |
| MeshTools::Generation                                                                                        |
|   build_cube()                   1         0.0113      0.011349    0.0113      0.011349    1.21     1.21     |
|                                                                                                              |
| MetisPartitioner                                                                                             |
|   partition()                    1         0.1345      0.134471    0.1721      0.172057    14.36    18.37    |
|                                                                                                              |
| Parallel                                                                                                     |
|   allgather()                    14        0.0018      0.000128    0.0020      0.000145    0.19     0.22     |
|   barrier()                      1         0.0016      0.001553    0.0016      0.001553    0.17     0.17     |
|   broadcast()                    13        0.0003      0.000025    0.0003      0.000021    0.03     0.03     |
|   max(bool)                      1         0.0000      0.000007    0.0000      0.000007    0.00     0.00     |
|   max(scalar)                    135       0.0021      0.000016    0.0021      0.000016    0.23     0.23     |
|   max(vector)                    30        0.0005      0.000018    0.0016      0.000053    0.06     0.17     |
|   min(bool)                      156       0.0017      0.000011    0.0017      0.000011    0.18     0.18     |
|   min(scalar)                    129       0.0267      0.000207    0.0267      0.000207    2.86     2.86     |
|   min(vector)                    30        0.0006      0.000022    0.0034      0.000112    0.07     0.36     |
|   probe()                        224       0.0050      0.000022    0.0050      0.000022    0.53     0.53     |
|   receive()                      202       0.0015      0.000008    0.0065      0.000032    0.16     0.69     |
|   send()                         202       0.0009      0.000005    0.0009      0.000005    0.10     0.10     |
|   send_receive()                 184       0.0015      0.000008    0.0088      0.000048    0.16     0.94     |
|   sum()                          31        0.0059      0.000192    0.0073      0.000235    0.63     0.78     |
|                                                                                                              |
| Parallel::Request                                                                                            |
|   wait()                         202       0.0004      0.000002    0.0004      0.000002    0.04     0.04     |
|                                                                                                              |
| Partitioner                                                                                                  |
|   set_node_processor_ids()       1         0.0048      0.004815    0.0081      0.008082    0.51     0.86     |
|   set_parent_processor_ids()     1         0.0044      0.004398    0.0044      0.004398    0.47     0.47     |
|                                                                                                              |
| RBConstruction                                                                                               |
|   clear()                        1         0.0007      0.000679    0.0007      0.000679    0.07     0.07     |
|   load_basis_function()          1         0.0001      0.000109    0.0001      0.000109    0.01     0.01     |
|   load_rb_solution()             1         0.0004      0.000421    0.0004      0.000421    0.04     0.04     |
|                                                                                                              |
| RBEvaluation                                                                                                 |
|   clear()                        1         0.0002      0.000202    0.0002      0.000202    0.02     0.02     |
|   compute_residual_dual_norm()   1         0.0025      0.002479    0.0025      0.002479    0.26     0.26     |
|   rb_solve()                     1         0.0113      0.011258    0.0137      0.013737    1.20     1.47     |
|   read_in_basis_functions()      1         0.0001      0.000078    0.1841      0.184080    0.01     19.66    |
|   read_in_vectors()              1         0.0756      0.075562    0.1840      0.184001    8.07     19.65    |
|   read_offline_data_from_files() 1         0.0017      0.001658    0.0024      0.002372    0.18     0.25     |
|   resize_data_structures()       1         0.0007      0.000714    0.0007      0.000714    0.08     0.08     |
 --------------------------------------------------------------------------------------------------------------
| Totals:                          4799      0.9366                                          100.00            |
 --------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example reduced_basis_ex1:
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
