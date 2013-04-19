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
</div>

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
          SimpleRBEvaluation(const Parallel::Communicator& comm)
            : TransientRBEvaluation(comm)
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
                                const std::string& name_in,
                                const unsigned int number_in)
          : Parent(es, name_in, number_in),
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
Build a mesh on the default MPI communicator.
</div>

<div class ="fragment">
<pre>
          Mesh mesh (init.comm(), dim);
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
          SimpleRBEvaluation rb_eval(mesh.comm());
        
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
    SimpleRBEvaluation(<B><FONT COLOR="#228B22">const</FONT></B> Parallel::Communicator&amp; comm)
      : TransientRBEvaluation(comm)
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
                          <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp; name_in,
                          <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> number_in)
    : Parent(es, name_in, number_in),
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
  
    Mesh mesh (init.comm(), dim);
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
  
    SimpleRBEvaluation rb_eval(mesh.comm());
  
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
make[4]: Entering directory `/net/spark/workspace/roystgnr/libmesh/git/devel/examples/reduced_basis/reduced_basis_ex3'
***************************************************************
* Running Example reduced_basis_ex3:
*  mpirun -np 4 example-devel -online_mode 0 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc
***************************************************************
 
*** Warning, This code is untested, experimental, or likely to see future API changes: ../src/reduced_basis/rb_parametrized.C, line 41, compiled Apr 19 2013 at 11:42:51 ***
 EquationSystems
  n_systems()=1
   System #0, "RBConvectionDiffusion"
    Type "TransientRBConstruction"
    Variables="u" 
    Finite Element Types="LAGRANGE", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" 
    Approximation Orders="FIRST", "THIRD" 
    n_dofs()=676
    n_local_dofs()=186
    n_constrained_dofs()=100
    n_local_constrained_dofs()=25
    n_vectors()=3
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 8.33136
      Average Off-Processor Bandwidth <= 0.488166
      Maximum  On-Processor Bandwidth <= 11
      Maximum Off-Processor Bandwidth <= 8
    DofMap Constraints
      Number of DoF Constraints = 100
      Average DoF Constraint Length= 0
      Number of Node Constraints = 0

 Mesh Information:
  mesh_dimension()=2
  spatial_dimension()=3
  n_nodes()=676
    n_local_nodes()=186
  n_elem()=625
    n_local_elem()=156
    n_active_elem()=625
  n_subdomains()=1
  n_partitions()=4
  n_processors()=4
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
last eigenvalue = -1.81536e-16

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
last eigenvalue = -3.81639e-17

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
last eigenvalue = -2.99661e-16

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
last eigenvalue = -5.39623e-17

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
last eigenvalue = -5.17985e-17

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
last eigenvalue = -1.00369e-17

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
last eigenvalue = -1.03771e-18

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
last eigenvalue = -4.31652e-17

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
last eigenvalue = -1.35639e-17

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
last eigenvalue = -8.529e-18

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
last eigenvalue = -2.97398e-17

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
last eigenvalue = -4.98733e-18

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
last eigenvalue = -6.83047e-18

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
last eigenvalue = -2.44655e-18

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
last eigenvalue = -7.08296e-19

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
last eigenvalue = -1.3397e-18

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
last eigenvalue = -2.38181e-19

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
last eigenvalue = -6.64074e-19

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
last eigenvalue = -3.99824e-19

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
last eigenvalue = -1.95156e-18

Updating RB matrices
Updating RB residual terms
Updating RB initial conditions

---- Basis dimension: 20 ----
Performing RB solves on training set
Maximum (relative) error bound is 0.145574

Maximum number of basis functions reached: Nmax = 20

 -------------------------------------------------------------------------------------------------------------------
| Processor id:   0                                                                                                 |
| Num Processors: 4                                                                                                 |
| Time:           Fri Apr 19 11:55:20 2013                                                                          |
| OS:             Linux                                                                                             |
| HostName:       spark.ices.utexas.edu                                                                             |
| OS Release:     2.6.32-279.22.1.el6.x86_64                                                                        |
| OS Version:     #1 SMP Tue Feb 5 14:33:39 CST 2013                                                                |
| Machine:        x86_64                                                                                            |
| Username:       roystgnr                                                                                          |
| Configuration:  ../configure  '--enable-everything'                                                               |
|  'METHODS=devel'                                                                                                  |
|  '--prefix=/h2/roystgnr/libmesh-test'                                                                             |
|  'CXX=distcc /usr/bin/g++'                                                                                        |
|  'CC=distcc /usr/bin/gcc'                                                                                         |
|  'FC=distcc /usr/bin/gfortran'                                                                                    |
|  'F77=distcc /usr/bin/gfortran'                                                                                   |
|  'PETSC_DIR=/opt/apps/ossw/libraries/petsc/petsc-3.3-p2'                                                          |
|  'PETSC_ARCH=gcc-system-mkl-gf-10.3.12.361-mpich2-1.4.1p1-cxx-opt'                                                |
|  'SLEPC_DIR=/opt/apps/ossw/libraries/slepc/slepc-3.3-p2-petsc-3.3-p2-cxx-opt'                                     |
|  'TRILINOS_DIR=/opt/apps/ossw/libraries/trilinos/trilinos-10.12.2/sl6/gcc-system/mpich2-1.4.1p1/mkl-gf-10.3.12.361'|
|  'VTK_DIR=/opt/apps/ossw/libraries/vtk/vtk-5.10.0/sl6/gcc-system'                                                 |
|  'HDF5_DIR=/opt/apps/ossw/libraries/hdf5/hdf5-1.8.9/sl6/gcc-system'                                               |
 -------------------------------------------------------------------------------------------------------------------
 -------------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=9.24779, Active time=9.04028                                                      |
 -------------------------------------------------------------------------------------------------------------------
| Event                                 nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                                 w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|-------------------------------------------------------------------------------------------------------------------|
|                                                                                                                   |
|                                                                                                                   |
| DofMap                                                                                                            |
|   add_neighbors_to_send_list()        1         0.0008      0.000782    0.0011      0.001114    0.01     0.01     |
|   build_constraint_matrix()           1716      0.0027      0.000002    0.0027      0.000002    0.03     0.03     |
|   build_sparsity()                    1         0.0007      0.000686    0.0018      0.001849    0.01     0.02     |
|   cnstrn_elem_mat_vec()               1716      0.0098      0.000006    0.0098      0.000006    0.11     0.11     |
|   create_dof_constraints()            1         0.0019      0.001903    0.0050      0.004989    0.02     0.06     |
|   distribute_dofs()                   1         0.0016      0.001552    0.0044      0.004396    0.02     0.05     |
|   dof_indices()                       4282      0.0209      0.000005    0.0209      0.000005    0.23     0.23     |
|   prepare_send_list()                 1         0.0000      0.000008    0.0000      0.000008    0.00     0.00     |
|   reinit()                            1         0.0025      0.002458    0.0025      0.002458    0.03     0.03     |
|                                                                                                                   |
| FE                                                                                                                |
|   compute_shape_functions()           3960      0.0180      0.000005    0.0180      0.000005    0.20     0.20     |
|   init_shape_functions()              550       0.0033      0.000006    0.0033      0.000006    0.04     0.04     |
|   inverse_map()                       1056      0.0037      0.000004    0.0037      0.000004    0.04     0.04     |
|                                                                                                                   |
| FEMap                                                                                                             |
|   compute_affine_map()                3960      0.0109      0.000003    0.0109      0.000003    0.12     0.12     |
|   compute_face_map()                  528       0.0033      0.000006    0.0071      0.000013    0.04     0.08     |
|   init_face_shape_functions()         22        0.0001      0.000003    0.0001      0.000003    0.00     0.00     |
|   init_reference_to_physical_map()    550       0.0025      0.000004    0.0025      0.000004    0.03     0.03     |
|                                                                                                                   |
| Mesh                                                                                                              |
|   find_neighbors()                    1         0.0019      0.001870    0.0020      0.001994    0.02     0.02     |
|   renumber_nodes_and_elem()           2         0.0002      0.000097    0.0002      0.000097    0.00     0.00     |
|                                                                                                                   |
| MeshCommunication                                                                                                 |
|   assign_global_indices()             1         0.0078      0.007798    0.0114      0.011413    0.09     0.13     |
|   compute_hilbert_indices()           2         0.0051      0.002555    0.0051      0.002555    0.06     0.06     |
|   find_global_indices()               2         0.0083      0.004138    0.0145      0.007234    0.09     0.16     |
|   parallel_sort()                     2         0.0006      0.000292    0.0007      0.000345    0.01     0.01     |
|                                                                                                                   |
| MeshTools::Generation                                                                                             |
|   build_cube()                        1         0.0017      0.001734    0.0017      0.001734    0.02     0.02     |
|                                                                                                                   |
| MetisPartitioner                                                                                                  |
|   partition()                         1         0.0052      0.005231    0.0160      0.016031    0.06     0.18     |
|                                                                                                                   |
| Parallel                                                                                                          |
|   allgather()                         14        0.0030      0.000213    0.0030      0.000216    0.03     0.03     |
|   barrier()                           1         0.0000      0.000018    0.0000      0.000018    0.00     0.00     |
|   broadcast()                         44        0.0002      0.000005    0.0002      0.000005    0.00     0.00     |
|   max(bool)                           7         0.0000      0.000005    0.0000      0.000005    0.00     0.00     |
|   max(scalar)                         127       0.0016      0.000012    0.0016      0.000012    0.02     0.02     |
|   max(vector)                         30        0.0004      0.000014    0.0014      0.000047    0.00     0.02     |
|   maxloc(scalar)                      21        0.3881      0.018479    0.3881      0.018479    4.29     4.29     |
|   min(bool)                           148       0.0017      0.000011    0.0017      0.000011    0.02     0.02     |
|   min(scalar)                         121       0.0039      0.000032    0.0039      0.000032    0.04     0.04     |
|   min(vector)                         30        0.0005      0.000017    0.0016      0.000055    0.01     0.02     |
|   probe()                             80        0.0015      0.000018    0.0015      0.000018    0.02     0.02     |
|   receive()                           64        0.0003      0.000004    0.0005      0.000008    0.00     0.01     |
|   send()                              52        0.0001      0.000002    0.0001      0.000002    0.00     0.00     |
|   send_receive()                      56        0.0003      0.000006    0.0009      0.000016    0.00     0.01     |
|   sum()                               44        0.0008      0.000018    0.0008      0.000019    0.01     0.01     |
|                                                                                                                   |
| Parallel::Request                                                                                                 |
|   wait()                              52        0.0001      0.000001    0.0001      0.000001    0.00     0.00     |
|                                                                                                                   |
| Partitioner                                                                                                       |
|   set_node_processor_ids()            1         0.0005      0.000537    0.0007      0.000724    0.01     0.01     |
|   set_parent_processor_ids()          1         0.0002      0.000181    0.0002      0.000181    0.00     0.00     |
|                                                                                                                   |
| PetscLinearSolver                                                                                                 |
|   solve()                             2085      1.3715      0.000658    1.3715      0.000658    15.17    15.17    |
|                                                                                                                   |
| RBConstruction                                                                                                    |
|   add_scaled_matrix_and_vector()      11        0.0336      0.003053    0.1064      0.009676    0.37     1.18     |
|   clear()                             3         0.0005      0.000167    0.0005      0.000167    0.01     0.01     |
|   compute_Fq_representor_innerprods() 1         0.0008      0.000761    0.0037      0.003687    0.01     0.04     |
|   compute_max_error_bound()           21        0.0058      0.000278    1.0921      0.052005    0.06     12.08    |
|   compute_output_dual_innerprods()    1         0.0026      0.002554    0.0085      0.008458    0.03     0.09     |
|   train_reduced_basis()               1         0.0036      0.003564    8.8675      8.867546    0.04     98.09    |
|   update_RB_system_matrices()         20        0.0548      0.002742    0.0548      0.002742    0.61     0.61     |
|   update_residual_terms()             20        0.1207      0.006034    0.2545      0.012723    1.33     2.81     |
|                                                                                                                   |
| RBEvaluation                                                                                                      |
|   clear()                             2         0.0004      0.000208    0.0004      0.000208    0.00     0.00     |
|   resize_data_structures()            1         0.0001      0.000097    0.0001      0.000097    0.00     0.00     |
|   write_offline_data_to_files()       1         0.0103      0.010260    0.0103      0.010260    0.11     0.11     |
|   write_out_basis_functions()         1         0.0001      0.000110    0.0191      0.019066    0.00     0.21     |
|   write_out_vectors()                 1         0.0052      0.005225    0.0190      0.018956    0.06     0.21     |
|                                                                                                                   |
| TransientRBConstruction                                                                                           |
|   enrich_RB_space()                   20        0.4092      0.020459    0.4092      0.020459    4.53     4.53     |
|   mass_matrix_scaled_matvec()         2000      0.2502      0.000125    0.2502      0.000125    2.77     2.77     |
|   set_error_temporal_data()           2020      0.3237      0.000160    0.3237      0.000160    3.58     3.58     |
|   truth_assembly()                    2000      4.8078      0.002404    5.0582      0.002529    53.18    55.95    |
|   truth_solve()                       20        0.3025      0.015125    6.8730      0.343649    3.35     76.03    |
|   update_RB_initial_condition_all_N() 20        0.0346      0.001732    0.0346      0.001732    0.38     0.38     |
|   update_RB_system_matrices()         20        0.0175      0.000873    0.0723      0.003615    0.19     0.80     |
|   update_residual_terms()             20        0.0749      0.003745    0.3705      0.018523    0.83     4.10     |
|                                                                                                                   |
| TransientRBEvaluation                                                                                             |
|   cache_online_residual_terms()       525       0.0086      0.000016    0.0086      0.000016    0.09     0.09     |
|   compute_residual_dual_norm()        52500     0.4414      0.000008    0.4414      0.000008    4.88     4.88     |
|   rb_solve()                          525       0.2431      0.000463    0.6977      0.001329    2.69     7.72     |
|   resize_data_structures()            1         0.0001      0.000069    0.0002      0.000166    0.00     0.00     |
|   write_offline_data_to_files()       1         0.0005      0.000538    0.0108      0.010798    0.01     0.12     |
 -------------------------------------------------------------------------------------------------------------------
| Totals:                               81092     9.0403                                          100.00            |
 -------------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example reduced_basis_ex3:
*  mpirun -np 4 example-devel -online_mode 0 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc
***************************************************************
***************************************************************
* Running Example reduced_basis_ex3:
*  mpirun -np 4 example-devel -online_mode 1 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc
***************************************************************
 
*** Warning, This code is untested, experimental, or likely to see future API changes: ../src/reduced_basis/rb_parametrized.C, line 41, compiled Apr 19 2013 at 11:42:51 ***
 EquationSystems
  n_systems()=1
   System #0, "RBConvectionDiffusion"
    Type "TransientRBConstruction"
    Variables="u" 
    Finite Element Types="LAGRANGE", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" 
    Approximation Orders="FIRST", "THIRD" 
    n_dofs()=676
    n_local_dofs()=186
    n_constrained_dofs()=100
    n_local_constrained_dofs()=25
    n_vectors()=3
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 8.33136
      Average Off-Processor Bandwidth <= 0.488166
      Maximum  On-Processor Bandwidth <= 11
      Maximum Off-Processor Bandwidth <= 8
    DofMap Constraints
      Number of DoF Constraints = 100
      Average DoF Constraint Length= 0
      Number of Node Constraints = 0

 Mesh Information:
  mesh_dimension()=2
  spatial_dimension()=3
  n_nodes()=676
    n_local_nodes()=186
  n_elem()=625
    n_local_elem()=156
    n_active_elem()=625
  n_subdomains()=1
  n_partitions()=4
  n_processors()=4
  n_threads()=1
  processor_id()=0

x_vel: 1.000000e+00
y_vel: 1.000000e+00

Error bound (absolute) at the final time is 0.0167619


 -------------------------------------------------------------------------------------------------------------------
| Processor id:   0                                                                                                 |
| Num Processors: 4                                                                                                 |
| Time:           Fri Apr 19 11:55:21 2013                                                                          |
| OS:             Linux                                                                                             |
| HostName:       spark.ices.utexas.edu                                                                             |
| OS Release:     2.6.32-279.22.1.el6.x86_64                                                                        |
| OS Version:     #1 SMP Tue Feb 5 14:33:39 CST 2013                                                                |
| Machine:        x86_64                                                                                            |
| Username:       roystgnr                                                                                          |
| Configuration:  ../configure  '--enable-everything'                                                               |
|  'METHODS=devel'                                                                                                  |
|  '--prefix=/h2/roystgnr/libmesh-test'                                                                             |
|  'CXX=distcc /usr/bin/g++'                                                                                        |
|  'CC=distcc /usr/bin/gcc'                                                                                         |
|  'FC=distcc /usr/bin/gfortran'                                                                                    |
|  'F77=distcc /usr/bin/gfortran'                                                                                   |
|  'PETSC_DIR=/opt/apps/ossw/libraries/petsc/petsc-3.3-p2'                                                          |
|  'PETSC_ARCH=gcc-system-mkl-gf-10.3.12.361-mpich2-1.4.1p1-cxx-opt'                                                |
|  'SLEPC_DIR=/opt/apps/ossw/libraries/slepc/slepc-3.3-p2-petsc-3.3-p2-cxx-opt'                                     |
|  'TRILINOS_DIR=/opt/apps/ossw/libraries/trilinos/trilinos-10.12.2/sl6/gcc-system/mpich2-1.4.1p1/mkl-gf-10.3.12.361'|
|  'VTK_DIR=/opt/apps/ossw/libraries/vtk/vtk-5.10.0/sl6/gcc-system'                                                 |
|  'HDF5_DIR=/opt/apps/ossw/libraries/hdf5/hdf5-1.8.9/sl6/gcc-system'                                               |
 -------------------------------------------------------------------------------------------------------------------
 --------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=0.226937, Active time=0.154002                                               |
 --------------------------------------------------------------------------------------------------------------
| Event                            nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                            w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|--------------------------------------------------------------------------------------------------------------|
|                                                                                                              |
|                                                                                                              |
| DofMap                                                                                                       |
|   add_neighbors_to_send_list()   1         0.0005      0.000467    0.0007      0.000664    0.30     0.43     |
|   build_sparsity()               1         0.0004      0.000401    0.0017      0.001652    0.26     1.07     |
|   create_dof_constraints()       1         0.0011      0.001084    0.0028      0.002844    0.70     1.85     |
|   distribute_dofs()              1         0.0009      0.000948    0.0038      0.003843    0.62     2.50     |
|   dof_indices()                  1162      0.0032      0.000003    0.0032      0.000003    2.08     2.08     |
|   prepare_send_list()            1         0.0000      0.000006    0.0000      0.000006    0.00     0.00     |
|   reinit()                       1         0.0015      0.001474    0.0015      0.001474    0.96     0.96     |
|                                                                                                              |
| EquationSystems                                                                                              |
|   build_solution_vector()        2         0.0006      0.000285    0.0045      0.002244    0.37     2.91     |
|                                                                                                              |
| ExodusII_IO                                                                                                  |
|   write_nodal_data()             2         0.0782      0.039091    0.0782      0.039091    50.77    50.77    |
|                                                                                                              |
| Mesh                                                                                                         |
|   find_neighbors()               1         0.0021      0.002103    0.0022      0.002194    1.37     1.42     |
|   renumber_nodes_and_elem()      2         0.0001      0.000075    0.0001      0.000075    0.10     0.10     |
|                                                                                                              |
| MeshCommunication                                                                                            |
|   assign_global_indices()        1         0.0136      0.013561    0.0140      0.013968    8.81     9.07     |
|   compute_hilbert_indices()      2         0.0040      0.002007    0.0040      0.002007    2.61     2.61     |
|   find_global_indices()          2         0.0005      0.000264    0.0062      0.003122    0.34     4.05     |
|   parallel_sort()                2         0.0004      0.000188    0.0013      0.000648    0.24     0.84     |
|                                                                                                              |
| MeshOutput                                                                                                   |
|   write_equation_systems()       2         0.0001      0.000038    0.0829      0.041451    0.05     53.83    |
|                                                                                                              |
| MeshTools::Generation                                                                                        |
|   build_cube()                   1         0.0006      0.000559    0.0006      0.000559    0.36     0.36     |
|                                                                                                              |
| MetisPartitioner                                                                                             |
|   partition()                    1         0.0035      0.003462    0.0067      0.006655    2.25     4.32     |
|                                                                                                              |
| Parallel                                                                                                     |
|   allgather()                    14        0.0010      0.000075    0.0011      0.000080    0.68     0.73     |
|   barrier()                      1         0.0000      0.000019    0.0000      0.000019    0.01     0.01     |
|   broadcast()                    13        0.0001      0.000005    0.0001      0.000004    0.05     0.04     |
|   max(bool)                      1         0.0000      0.000002    0.0000      0.000002    0.00     0.00     |
|   max(scalar)                    135       0.0023      0.000017    0.0023      0.000017    1.46     1.46     |
|   max(vector)                    30        0.0002      0.000008    0.0007      0.000024    0.15     0.47     |
|   min(bool)                      156       0.0008      0.000005    0.0008      0.000005    0.53     0.53     |
|   min(scalar)                    129       0.0075      0.000058    0.0075      0.000058    4.87     4.87     |
|   min(vector)                    30        0.0003      0.000010    0.0009      0.000030    0.19     0.58     |
|   probe()                        64        0.0006      0.000010    0.0006      0.000010    0.40     0.40     |
|   receive()                      58        0.0002      0.000003    0.0008      0.000013    0.11     0.51     |
|   send()                         58        0.0002      0.000003    0.0002      0.000003    0.11     0.11     |
|   send_receive()                 56        0.0002      0.000004    0.0011      0.000019    0.15     0.71     |
|   sum()                          31        0.0011      0.000036    0.0041      0.000132    0.72     2.66     |
|                                                                                                              |
| Parallel::Request                                                                                            |
|   wait()                         58        0.0001      0.000001    0.0001      0.000001    0.04     0.04     |
|                                                                                                              |
| Partitioner                                                                                                  |
|   set_node_processor_ids()       1         0.0003      0.000320    0.0018      0.001812    0.21     1.18     |
|   set_parent_processor_ids()     1         0.0001      0.000103    0.0001      0.000103    0.07     0.07     |
|                                                                                                              |
| RBConstruction                                                                                               |
|   clear()                        3         0.0002      0.000068    0.0002      0.000068    0.13     0.13     |
|   load_basis_function()          1         0.0001      0.000106    0.0001      0.000106    0.07     0.07     |
|                                                                                                              |
| RBEvaluation                                                                                                 |
|   clear()                        2         0.0000      0.000018    0.0000      0.000018    0.02     0.02     |
|   read_in_basis_functions()      1         0.0000      0.000015    0.0208      0.020762    0.01     13.48    |
|   read_in_vectors()              1         0.0061      0.006135    0.0207      0.020747    3.98     13.47    |
|   read_offline_data_from_files() 1         0.0007      0.000688    0.0008      0.000815    0.45     0.53     |
|   resize_data_structures()       1         0.0001      0.000088    0.0001      0.000088    0.06     0.06     |
|                                                                                                              |
| TransientRBConstruction                                                                                      |
|   load_rb_solution()             1         0.0001      0.000146    0.0001      0.000146    0.09     0.09     |
|                                                                                                              |
| TransientRBEvaluation                                                                                        |
|   cache_online_residual_terms()  1         0.0000      0.000043    0.0000      0.000043    0.03     0.03     |
|   compute_residual_dual_norm()   100       0.0039      0.000039    0.0039      0.000039    2.54     2.54     |
|   rb_solve()                     1         0.0161      0.016114    0.0201      0.020084    10.46    13.04    |
|   read_offline_data_from_files() 1         0.0003      0.000289    0.0011      0.001104    0.19     0.72     |
|   resize_data_structures()       1         0.0000      0.000039    0.0001      0.000127    0.03     0.08     |
 --------------------------------------------------------------------------------------------------------------
| Totals:                          2137      0.1540                                          100.00            |
 --------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example reduced_basis_ex3:
*  mpirun -np 4 example-devel -online_mode 1 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc
***************************************************************
make[4]: Leaving directory `/net/spark/workspace/roystgnr/libmesh/git/devel/examples/reduced_basis/reduced_basis_ex3'
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
