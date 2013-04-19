<?php $root=""; ?>
<?php require($root."navigation.php"); ?>
<html>
<head>
  <?php load_style($root); ?>
</head>
 
<body>
 
<?php make_navigation("reduced_basis_ex4",$root)?>
 
<div class="content">
<a name="comments"></a> 
<br><br><br> <h1> The source file assembly.h with comments: </h1> 
<div class = "comment">
</div>

<div class ="fragment">
<pre>
        #ifndef __assembly_h__
        #define __assembly_h__
        
</pre>
</div>
<div class = "comment">
libMesh includes
</div>

<div class ="fragment">
<pre>
        #include "libmesh/libmesh.h"
        #include "libmesh/mesh.h"
        #include "libmesh/equation_systems.h"
        #include "libmesh/fe.h"
        #include "libmesh/quadrature.h"
        #include "libmesh/dof_map.h"
        #include "libmesh/dense_matrix.h"
        #include "libmesh/dense_vector.h"
        #include "libmesh/fe_interface.h"
        #include "libmesh/elem.h"
        
</pre>
</div>
<div class = "comment">
rbOOmit includes
</div>

<div class ="fragment">
<pre>
        #include "libmesh/rb_assembly_expansion.h"
        #include "libmesh/rb_eim_theta.h"
        #include "libmesh/rb_parametrized_function.h"
        
</pre>
</div>
<div class = "comment">
Bring in bits from the libMesh namespace.
Just the bits we're using, since this is a header.
</div>

<div class ="fragment">
<pre>
        using libMesh::ElemAssembly;
        using libMesh::FEMContext;
        using libMesh::Number;
        using libMesh::Point;
        using libMesh::RBAssemblyExpansion;
        using libMesh::RBEIMAssembly;
        using libMesh::RBEIMConstruction;
        using libMesh::RBParametrizedFunction;
        using libMesh::RBParameters;
        using libMesh::RBTheta;
        using libMesh::RBThetaExpansion;
        using libMesh::Real;
        using libMesh::RealGradient;
        
        struct ShiftedGaussian : public RBParametrizedFunction
        {
          virtual Number evaluate(const RBParameters& mu,
                                  const Point& p)
          {
            Real center_x = mu.get_value("center_x");
            Real center_y = mu.get_value("center_y");
            return exp( -2.*(pow(center_x-p(0),2.) + pow(center_y-p(1),2.)) );
          }
        };
        
</pre>
</div>
<div class = "comment">
Expansion of the PDE operator
</div>

<div class ="fragment">
<pre>
        struct ThetaA0 : RBTheta { virtual Number evaluate(const RBParameters& )   { return 0.05;  } };
        
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
        
        
        struct EIM_IP_assembly : ElemAssembly
        {
        
</pre>
</div>
<div class = "comment">
Use the L2 norm to find the best fit
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
        
            const unsigned int n_u_dofs = c.dof_indices_var[u_var].size();
        
            unsigned int n_qpoints = (c.get_element_qrule())-&gt;n_points();
        
            for (unsigned int qp=0; qp != n_qpoints; qp++)
              for (unsigned int i=0; i != n_u_dofs; i++)
                for (unsigned int j=0; j != n_u_dofs; j++)
                  c.get_elem_jacobian()(i,j) += JxW[qp] * phi[j][qp]*phi[i][qp];
          }
        };
        
        struct EIM_F : RBEIMAssembly
        {
        
          EIM_F(RBEIMConstruction& rb_eim_con_in,
                unsigned int basis_function_index_in)
          : RBEIMAssembly(rb_eim_con_in,
                          basis_function_index_in)
          {}
        
          virtual void interior_assembly(FEMContext &c)
          {
</pre>
</div>
<div class = "comment">
PDE variable number
</div>

<div class ="fragment">
<pre>
            const unsigned int u_var = 0;
        
</pre>
</div>
<div class = "comment">
EIM variable number
</div>

<div class ="fragment">
<pre>
            const unsigned int eim_var = 0;
        
            const std::vector&lt;Real&gt; &JxW =
              c.element_fe_var[u_var]-&gt;get_JxW();
        
            const std::vector&lt;std::vector&lt;Real&gt; &gt;& phi =
              c.element_fe_var[u_var]-&gt;get_phi();
        
            const std::vector&lt;Point&gt;& qpoints =
              c.element_fe_var[u_var]-&gt;get_xyz();
        
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
        
            std::vector&lt;Number&gt; eim_values;
            evaluate_basis_function(eim_var,
                                    *c.elem,
                                    qpoints,
                                    eim_values);
        
            for (unsigned int qp=0; qp != n_qpoints; qp++)
              for (unsigned int i=0; i != n_u_dofs; i++)
                c.get_elem_residual()(i) += JxW[qp] * ( eim_values[qp]*phi[i][qp] );
          }
        
        };
        
</pre>
</div>
<div class = "comment">
Define an RBThetaExpansion class for this PDE
</div>

<div class ="fragment">
<pre>
        struct EimTestRBThetaExpansion : RBThetaExpansion
        {
        
          /**
           * Constructor.
           */
          EimTestRBThetaExpansion()
          {
            attach_A_theta(&theta_a_0);
          }
        
</pre>
</div>
<div class = "comment">
The RBTheta member variables
</div>

<div class ="fragment">
<pre>
          ThetaA0 theta_a_0;
        };
        
</pre>
</div>
<div class = "comment">
Define an RBAssemblyExpansion class for this PDE
</div>

<div class ="fragment">
<pre>
        struct EimTestRBAssemblyExpansion : RBAssemblyExpansion
        {
        
          /**
           * Constructor.
           */
          EimTestRBAssemblyExpansion()
          {
            attach_A_assembly(&A0_assembly);
          }
        
</pre>
</div>
<div class = "comment">
A0 assembly object
</div>

<div class ="fragment">
<pre>
          A0 A0_assembly;
        
        };
        
        #endif
</pre>
</div>

<a name="comments"></a> 
<br><br><br> <h1> The source file eim_classes.h with comments: </h1> 
<div class = "comment">
</div>

<div class ="fragment">
<pre>
        #ifndef __eim_classes_h__
        #define __eim_classes_h__
        
</pre>
</div>
<div class = "comment">
local includes
</div>

<div class ="fragment">
<pre>
        #include "libmesh/rb_eim_construction.h"
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
        using libMesh::EquationSystems;
        using libMesh::RBEIMEvaluation;
        
</pre>
</div>
<div class = "comment">
A simple subclass of RBEIMEvaluation. Overload
evaluate_parametrized_function to define the
function that we "empirically" interpolate.
</div>

<div class ="fragment">
<pre>
        class SimpleEIMEvaluation : public RBEIMEvaluation
        {
        public:
        
          SimpleEIMEvaluation(const Parallel::Communicator& comm)
            : RBEIMEvaluation(comm)
          {
            attach_parametrized_function(&sg);
          }
        
          /**
           * Parametrized function that we approximate with EIM
           */
          ShiftedGaussian sg;
        
        };
        
</pre>
</div>
<div class = "comment">
A simple subclass of RBEIMConstruction.
</div>

<div class ="fragment">
<pre>
        class SimpleEIMConstruction : public RBEIMConstruction
        {
        public:
        
          /**
           * Constructor.
           */
          SimpleEIMConstruction (EquationSystems& es,
                                 const std::string& name_in,
                                 const unsigned int number_in)
          : Parent(es, name_in, number_in)
          {
          }
        
          /**
           * The type of the parent.
           */
          typedef RBEIMConstruction Parent;
        
          /**
           * Provide an implementation of build_eim_assembly
           */
          virtual AutoPtr&lt;ElemAssembly&gt; build_eim_assembly(unsigned int index)
          {
            return AutoPtr&lt;ElemAssembly&gt;(new EIM_F(*this, index));
          }
        
          /**
           * Initialize data structures.
           */
          virtual void init_data()
          {
            u_var = this-&gt;add_variable ("f_EIM", FIRST);
        
            Parent::init_data();
        
            set_inner_product_assembly(ip);
          }
        
          /**
           * Variable number for u.
           */
          unsigned int u_var;
        
          /**
           * Inner product assembly object
           */
          EIM_IP_assembly ip;
        
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
        
        #include "libmesh/rb_construction.h"
        #include "assembly.h"
        
</pre>
</div>
<div class = "comment">
Bring in bits from the libMesh namespace.
Just the bits we're using, since this is a header.
</div>

<div class ="fragment">
<pre>
        using libMesh::DirichletBoundary;
        using libMesh::RBConstruction;
        using libMesh::RBEvaluation;
        
</pre>
</div>
<div class = "comment">
A simple subclass of RBEvaluation.
</div>

<div class ="fragment">
<pre>
        class SimpleRBEvaluation : public RBEvaluation
        {
        public:
        
          /**
           * Constructor. Just set the theta expansion.
           */
          SimpleRBEvaluation(const Parallel::Communicator& comm)
            : RBEvaluation(comm)
          {
            set_rb_theta_expansion(eim_test_rb_theta_expansion);
          }
        
          /**
           * The object that stores the "theta" expansion of the parameter dependent PDE,
           * i.e. the set of parameter-dependent functions in the affine expansion of the PDE.
           */
          EimTestRBThetaExpansion eim_test_rb_theta_expansion;
        
        };
        
        class SimpleRBConstruction : public RBConstruction
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
</div>

<div class ="fragment">
<pre>
            set_rb_assembly_expansion(eim_test_rb_assembly_expansion);
        
</pre>
</div>
<div class = "comment">
We need to define an inner product matrix for this problem
</div>

<div class ="fragment">
<pre>
            set_inner_product_assembly(eim_test_rb_assembly_expansion.A0_assembly);
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
           * The object that stores the "theta" expansion of the parameter dependent PDE,
           * i.e. the set of parameter-dependent functions in the affine expansion of the PDE.
           */
          EimTestRBThetaExpansion eim_test_rb_theta_expansion;
        
          /**
           * The object that stores the "assembly" expansion of the parameter dependent PDE,
           * i.e. the objects that define how to assemble the set of parameter-independent
           * operators in the affine expansion of the PDE.
           */
          EimTestRBAssemblyExpansion eim_test_rb_assembly_expansion;
        
          /**
           * The object that defines which degrees of freedom are on a Dirichlet boundary.
           */
          AutoPtr&lt;DirichletBoundary&gt; dirichlet_bc;
        
        };
        
        #endif
</pre>
</div>

<a name="comments"></a> 
<br><br><br> <h1> The source file reduced_basis_ex4.C with comments: </h1> 
<div class = "comment">
</div>

<div class ="fragment">
<pre>
        #include "libmesh/libmesh.h"
        #include "libmesh/mesh.h"
        #include "libmesh/mesh_generation.h"
        #include "libmesh/equation_systems.h"
        #include "libmesh/exodusII_io.h"
        #include "libmesh/getpot.h"
        
        #include "eim_classes.h"
        #include "rb_classes.h"
        
        using namespace libMesh;
        
        
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
Skip this 2D example if libMesh was compiled as 1D-only.
</div>

<div class ="fragment">
<pre>
          libmesh_example_assert(2 &lt;= LIBMESH_DIM, "2D support");
        
</pre>
</div>
<div class = "comment">
Define the names of the input files we will read the problem properties from
</div>

<div class ="fragment">
<pre>
          std::string eim_parameters = "eim.in";
          std::string rb_parameters  = "rb.in";
          std::string main_parameters = "reduced_basis_ex4.in";
          GetPot infile(main_parameters);
        
          unsigned int n_elem = infile("n_elem", 1);       // Determines the number of elements in the "truth" mesh
          const unsigned int dim = 2;                      // The number of spatial dimensions
          bool store_basis_functions = infile("store_basis_functions", false); // Do we write out basis functions?
        
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
Create a mesh (just a simple square) on the default MPI
communicator
</div>

<div class ="fragment">
<pre>
          Mesh mesh (init.comm(), dim);
          MeshTools::Generation::build_square (mesh,
                                               n_elem, n_elem,
                                               -1., 1.,
                                               -1., 1.,
                                               QUAD4);
        
</pre>
</div>
<div class = "comment">
Initialize the EquationSystems object for this mesh and attach
the EIM and RB Construction objects
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
          mesh.print_info();
          equation_systems.print_info();
        
</pre>
</div>
<div class = "comment">
Initialize the standard RBEvaluation object
</div>

<div class ="fragment">
<pre>
          SimpleRBEvaluation rb_eval(mesh.comm());
        
</pre>
</div>
<div class = "comment">
Initialize the EIM RBEvaluation object
</div>

<div class ="fragment">
<pre>
          SimpleEIMEvaluation eim_rb_eval(mesh.comm());
        
</pre>
</div>
<div class = "comment">
Set the rb_eval objects for the RBConstructions
</div>

<div class ="fragment">
<pre>
          eim_construction.set_rb_evaluation(eim_rb_eval);
          rb_construction.set_rb_evaluation(rb_eval);
        
          if(!online_mode)
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
attach the EIM theta objects to the RBConstruction and RBEvaluation objects
</div>

<div class ="fragment">
<pre>
            eim_rb_eval.initialize_eim_theta_objects();
            rb_eval.get_rb_theta_expansion().attach_multiple_F_theta(eim_rb_eval.get_eim_theta_objects());
        
</pre>
</div>
<div class = "comment">
attach the EIM assembly objects to the RBConstruction object
</div>

<div class ="fragment">
<pre>
            eim_construction.initialize_eim_assembly_objects();
            rb_construction.get_rb_assembly_expansion().attach_multiple_F_assembly(eim_construction.get_eim_assembly_objects());
        
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
          else
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
            rb_eval.get_rb_theta_expansion().attach_multiple_F_theta(eim_rb_eval.get_eim_theta_objects());
        
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
            Real online_center_x = infile("online_center_x", 0.);
            Real online_center_y = infile("online_center_y", 0.);
            RBParameters online_mu;
            online_mu.set_value("center_x", online_center_x);
            online_mu.set_value("center_y", online_center_y);
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
        #ifdef LIBMESH_HAVE_EXODUS_API
              ExodusII_IO(mesh).write_equation_systems("RB_sol.e",equation_systems);
        #endif
            }
          }
        
        }
</pre>
</div>

<a name="nocomments"></a> 
<br><br><br> <h1> The source file assembly.h without comments: </h1> 
<pre> 
  #ifndef __assembly_h__
  #define __assembly_h__
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/libmesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/mesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/equation_systems.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/fe.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/quadrature.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/dof_map.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/dense_matrix.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/dense_vector.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/fe_interface.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/elem.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/rb_assembly_expansion.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/rb_eim_theta.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/rb_parametrized_function.h&quot;</FONT></B>
  
  using libMesh::ElemAssembly;
  using libMesh::FEMContext;
  using libMesh::Number;
  using libMesh::Point;
  using libMesh::RBAssemblyExpansion;
  using libMesh::RBEIMAssembly;
  using libMesh::RBEIMConstruction;
  using libMesh::RBParametrizedFunction;
  using libMesh::RBParameters;
  using libMesh::RBTheta;
  using libMesh::RBThetaExpansion;
  using libMesh::Real;
  using libMesh::RealGradient;
  
  <B><FONT COLOR="#228B22">struct</FONT></B> ShiftedGaussian : <B><FONT COLOR="#228B22">public</FONT></B> RBParametrizedFunction
  {
    <B><FONT COLOR="#228B22">virtual</FONT></B> Number evaluate(<B><FONT COLOR="#228B22">const</FONT></B> RBParameters&amp; mu,
                            <B><FONT COLOR="#228B22">const</FONT></B> Point&amp; p)
    {
      Real center_x = mu.get_value(<B><FONT COLOR="#BC8F8F">&quot;center_x&quot;</FONT></B>);
      Real center_y = mu.get_value(<B><FONT COLOR="#BC8F8F">&quot;center_y&quot;</FONT></B>);
      <B><FONT COLOR="#A020F0">return</FONT></B> exp( -2.*(pow(center_x-p(0),2.) + pow(center_y-p(1),2.)) );
    }
  };
  
  <B><FONT COLOR="#228B22">struct</FONT></B> ThetaA0 : RBTheta { <B><FONT COLOR="#228B22">virtual</FONT></B> Number evaluate(<B><FONT COLOR="#228B22">const</FONT></B> RBParameters&amp; )   { <B><FONT COLOR="#A020F0">return</FONT></B> 0.05;  } };
  
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
  
  
  <B><FONT COLOR="#228B22">struct</FONT></B> EIM_IP_assembly : ElemAssembly
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
            c.get_elem_jacobian()(i,j) += JxW[qp] * phi[j][qp]*phi[i][qp];
    }
  };
  
  <B><FONT COLOR="#228B22">struct</FONT></B> EIM_F : RBEIMAssembly
  {
  
    EIM_F(RBEIMConstruction&amp; rb_eim_con_in,
          <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> basis_function_index_in)
    : RBEIMAssembly(rb_eim_con_in,
                    basis_function_index_in)
    {}
  
    <B><FONT COLOR="#228B22">virtual</FONT></B> <B><FONT COLOR="#228B22">void</FONT></B> interior_assembly(FEMContext &amp;c)
    {
      <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> u_var = 0;
  
      <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> eim_var = 0;
  
      <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Real&gt; &amp;JxW =
        c.element_fe_var[u_var]-&gt;get_JxW();
  
      <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;Real&gt; &gt;&amp; phi =
        c.element_fe_var[u_var]-&gt;get_phi();
  
      <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Point&gt;&amp; qpoints =
        c.element_fe_var[u_var]-&gt;get_xyz();
  
      <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_u_dofs = c.dof_indices_var[u_var].size();
  
      <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_qpoints = (c.get_element_qrule())-&gt;n_points();
  
      <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;Number&gt; eim_values;
      evaluate_basis_function(eim_var,
                              *c.elem,
                              qpoints,
                              eim_values);
  
      <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> qp=0; qp != n_qpoints; qp++)
        <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i != n_u_dofs; i++)
          c.get_elem_residual()(i) += JxW[qp] * ( eim_values[qp]*phi[i][qp] );
    }
  
  };
  
  <B><FONT COLOR="#228B22">struct</FONT></B> EimTestRBThetaExpansion : RBThetaExpansion
  {
  
    <I><FONT COLOR="#B22222">/**
     * Constructor.
     */</FONT></I>
    EimTestRBThetaExpansion()
    {
      attach_A_theta(&amp;theta_a_0);
    }
  
    ThetaA0 theta_a_0;
  };
  
  <B><FONT COLOR="#228B22">struct</FONT></B> EimTestRBAssemblyExpansion : RBAssemblyExpansion
  {
  
    <I><FONT COLOR="#B22222">/**
     * Constructor.
     */</FONT></I>
    EimTestRBAssemblyExpansion()
    {
      attach_A_assembly(&amp;A0_assembly);
    }
  
    A0 A0_assembly;
  
  };
  
  #endif
</pre> 
<a name="nocomments"></a> 
<br><br><br> <h1> The source file eim_classes.h without comments: </h1> 
<pre> 
  #ifndef __eim_classes_h__
  #define __eim_classes_h__
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/rb_eim_construction.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;assembly.h&quot;</FONT></B>
  
  using libMesh::AutoPtr;
  using libMesh::EquationSystems;
  using libMesh::RBEIMEvaluation;
  
  <B><FONT COLOR="#228B22">class</FONT></B> SimpleEIMEvaluation : <B><FONT COLOR="#228B22">public</FONT></B> RBEIMEvaluation
  {
  <B><FONT COLOR="#228B22">public</FONT></B>:
  
    SimpleEIMEvaluation(<B><FONT COLOR="#228B22">const</FONT></B> Parallel::Communicator&amp; comm)
      : RBEIMEvaluation(comm)
    {
      attach_parametrized_function(&amp;sg);
    }
  
    <I><FONT COLOR="#B22222">/**
     * Parametrized function that we approximate with EIM
     */</FONT></I>
    ShiftedGaussian sg;
  
  };
  
  <B><FONT COLOR="#228B22">class</FONT></B> SimpleEIMConstruction : <B><FONT COLOR="#228B22">public</FONT></B> RBEIMConstruction
  {
  <B><FONT COLOR="#228B22">public</FONT></B>:
  
    <I><FONT COLOR="#B22222">/**
     * Constructor.
     */</FONT></I>
    SimpleEIMConstruction (EquationSystems&amp; es,
                           <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp; name_in,
                           <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> number_in)
    : Parent(es, name_in, number_in)
    {
    }
  
    <I><FONT COLOR="#B22222">/**
     * The type of the parent.
     */</FONT></I>
    <B><FONT COLOR="#228B22">typedef</FONT></B> RBEIMConstruction Parent;
  
    <I><FONT COLOR="#B22222">/**
     * Provide an implementation of build_eim_assembly
     */</FONT></I>
    <B><FONT COLOR="#228B22">virtual</FONT></B> AutoPtr&lt;ElemAssembly&gt; build_eim_assembly(<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> index)
    {
      <B><FONT COLOR="#A020F0">return</FONT></B> AutoPtr&lt;ElemAssembly&gt;(<B><FONT COLOR="#A020F0">new</FONT></B> EIM_F(*<B><FONT COLOR="#A020F0">this</FONT></B>, index));
    }
  
    <I><FONT COLOR="#B22222">/**
     * Initialize data structures.
     */</FONT></I>
    <B><FONT COLOR="#228B22">virtual</FONT></B> <B><FONT COLOR="#228B22">void</FONT></B> init_data()
    {
      u_var = <B><FONT COLOR="#A020F0">this</FONT></B>-&gt;add_variable (<B><FONT COLOR="#BC8F8F">&quot;f_EIM&quot;</FONT></B>, FIRST);
  
      <B><FONT COLOR="#5F9EA0">Parent</FONT></B>::init_data();
  
      set_inner_product_assembly(ip);
    }
  
    <I><FONT COLOR="#B22222">/**
     * Variable number for u.
     */</FONT></I>
    <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> u_var;
  
    <I><FONT COLOR="#B22222">/**
     * Inner product assembly object
     */</FONT></I>
    EIM_IP_assembly ip;
  
  };
  
  #endif
</pre> 
<a name="nocomments"></a> 
<br><br><br> <h1> The source file rb_classes.h without comments: </h1> 
<pre> 
  #ifndef __rb_classes_h__
  #define __rb_classes_h__
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/rb_construction.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;assembly.h&quot;</FONT></B>
  
  using libMesh::DirichletBoundary;
  using libMesh::RBConstruction;
  using libMesh::RBEvaluation;
  
  <B><FONT COLOR="#228B22">class</FONT></B> SimpleRBEvaluation : <B><FONT COLOR="#228B22">public</FONT></B> RBEvaluation
  {
  <B><FONT COLOR="#228B22">public</FONT></B>:
  
    <I><FONT COLOR="#B22222">/**
     * Constructor. Just set the theta expansion.
     */</FONT></I>
    SimpleRBEvaluation(<B><FONT COLOR="#228B22">const</FONT></B> Parallel::Communicator&amp; comm)
      : RBEvaluation(comm)
    {
      set_rb_theta_expansion(eim_test_rb_theta_expansion);
    }
  
    <I><FONT COLOR="#B22222">/**
     * The object that stores the &quot;theta&quot; expansion of the parameter dependent PDE,
     * i.e. the set of parameter-dependent functions in the affine expansion of the PDE.
     */</FONT></I>
    EimTestRBThetaExpansion eim_test_rb_theta_expansion;
  
  };
  
  <B><FONT COLOR="#228B22">class</FONT></B> SimpleRBConstruction : <B><FONT COLOR="#228B22">public</FONT></B> RBConstruction
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
  
      set_rb_assembly_expansion(eim_test_rb_assembly_expansion);
  
      set_inner_product_assembly(eim_test_rb_assembly_expansion.A0_assembly);
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
     * The object that stores the &quot;theta&quot; expansion of the parameter dependent PDE,
     * i.e. the set of parameter-dependent functions in the affine expansion of the PDE.
     */</FONT></I>
    EimTestRBThetaExpansion eim_test_rb_theta_expansion;
  
    <I><FONT COLOR="#B22222">/**
     * The object that stores the &quot;assembly&quot; expansion of the parameter dependent PDE,
     * i.e. the objects that define how to assemble the set of parameter-independent
     * operators in the affine expansion of the PDE.
     */</FONT></I>
    EimTestRBAssemblyExpansion eim_test_rb_assembly_expansion;
  
    <I><FONT COLOR="#B22222">/**
     * The object that defines which degrees of freedom are on a Dirichlet boundary.
     */</FONT></I>
    AutoPtr&lt;DirichletBoundary&gt; dirichlet_bc;
  
  };
  
  #endif
</pre> 
<a name="nocomments"></a> 
<br><br><br> <h1> The source file reduced_basis_ex4.C without comments: </h1> 
<pre> 
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/libmesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/mesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/mesh_generation.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/equation_systems.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/exodusII_io.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/getpot.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;eim_classes.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;rb_classes.h&quot;</FONT></B>
  
  using namespace libMesh;
  
  
  <B><FONT COLOR="#228B22">int</FONT></B> main (<B><FONT COLOR="#228B22">int</FONT></B> argc, <B><FONT COLOR="#228B22">char</FONT></B>** argv)
  {
    LibMeshInit init (argc, argv);
  
  #<B><FONT COLOR="#A020F0">if</FONT></B> !defined(LIBMESH_HAVE_XDR)
    libmesh_example_assert(false, <B><FONT COLOR="#BC8F8F">&quot;--enable-xdr&quot;</FONT></B>);
  #elif defined(LIBMESH_DEFAULT_SINGLE_PRECISION)
    libmesh_example_assert(false, <B><FONT COLOR="#BC8F8F">&quot;--disable-singleprecision&quot;</FONT></B>);
  #endif
  
    libmesh_example_assert(2 &lt;= LIBMESH_DIM, <B><FONT COLOR="#BC8F8F">&quot;2D support&quot;</FONT></B>);
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::string eim_parameters = <B><FONT COLOR="#BC8F8F">&quot;eim.in&quot;</FONT></B>;
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::string rb_parameters  = <B><FONT COLOR="#BC8F8F">&quot;rb.in&quot;</FONT></B>;
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::string main_parameters = <B><FONT COLOR="#BC8F8F">&quot;reduced_basis_ex4.in&quot;</FONT></B>;
    GetPot infile(main_parameters);
  
    <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_elem = infile(<B><FONT COLOR="#BC8F8F">&quot;n_elem&quot;</FONT></B>, 1);       <I><FONT COLOR="#B22222">// Determines the number of elements in the &quot;truth&quot; mesh
</FONT></I>    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> dim = 2;                      <I><FONT COLOR="#B22222">// The number of spatial dimensions
</FONT></I>    <B><FONT COLOR="#228B22">bool</FONT></B> store_basis_functions = infile(<B><FONT COLOR="#BC8F8F">&quot;store_basis_functions&quot;</FONT></B>, false); <I><FONT COLOR="#B22222">// Do we write out basis functions?
</FONT></I>  
    GetPot command_line (argc, argv);
    <B><FONT COLOR="#228B22">int</FONT></B> online_mode = 0;
    <B><FONT COLOR="#A020F0">if</FONT></B> ( command_line.search(1, <B><FONT COLOR="#BC8F8F">&quot;-online_mode&quot;</FONT></B>) )
      online_mode = command_line.next(online_mode);
  
    Mesh mesh (init.comm(), dim);
    <B><FONT COLOR="#5F9EA0">MeshTools</FONT></B>::Generation::build_square (mesh,
                                         n_elem, n_elem,
                                         -1., 1.,
                                         -1., 1.,
                                         QUAD4);
  
    EquationSystems equation_systems (mesh);
  
    SimpleEIMConstruction &amp; eim_construction =
      equation_systems.add_system&lt;SimpleEIMConstruction&gt; (<B><FONT COLOR="#BC8F8F">&quot;EIM&quot;</FONT></B>);
    SimpleRBConstruction &amp; rb_construction =
      equation_systems.add_system&lt;SimpleRBConstruction&gt; (<B><FONT COLOR="#BC8F8F">&quot;RB&quot;</FONT></B>);
  
    equation_systems.init ();
  
    mesh.print_info();
    equation_systems.print_info();
  
    SimpleRBEvaluation rb_eval(mesh.comm());
  
    SimpleEIMEvaluation eim_rb_eval(mesh.comm());
  
    eim_construction.set_rb_evaluation(eim_rb_eval);
    rb_construction.set_rb_evaluation(rb_eval);
  
    <B><FONT COLOR="#A020F0">if</FONT></B>(!online_mode)
    {
      eim_construction.process_parameters_file(eim_parameters);
      eim_construction.print_info();
  
      eim_construction.initialize_rb_construction();
      eim_construction.train_reduced_basis();
      eim_construction.get_rb_evaluation().write_offline_data_to_files(<B><FONT COLOR="#BC8F8F">&quot;eim_data&quot;</FONT></B>);
  
      rb_construction.process_parameters_file(rb_parameters);
  
      eim_rb_eval.initialize_eim_theta_objects();
      rb_eval.get_rb_theta_expansion().attach_multiple_F_theta(eim_rb_eval.get_eim_theta_objects());
  
      eim_construction.initialize_eim_assembly_objects();
      rb_construction.get_rb_assembly_expansion().attach_multiple_F_assembly(eim_construction.get_eim_assembly_objects());
  
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
    <B><FONT COLOR="#A020F0">else</FONT></B>
    {
      eim_rb_eval.read_offline_data_from_files(<B><FONT COLOR="#BC8F8F">&quot;eim_data&quot;</FONT></B>);
  
      eim_rb_eval.initialize_eim_theta_objects();
      rb_eval.get_rb_theta_expansion().attach_multiple_F_theta(eim_rb_eval.get_eim_theta_objects());
  
      rb_eval.read_offline_data_from_files(<B><FONT COLOR="#BC8F8F">&quot;rb_data&quot;</FONT></B>);
  
      Real online_center_x = infile(<B><FONT COLOR="#BC8F8F">&quot;online_center_x&quot;</FONT></B>, 0.);
      Real online_center_y = infile(<B><FONT COLOR="#BC8F8F">&quot;online_center_y&quot;</FONT></B>, 0.);
      RBParameters online_mu;
      online_mu.set_value(<B><FONT COLOR="#BC8F8F">&quot;center_x&quot;</FONT></B>, online_center_x);
      online_mu.set_value(<B><FONT COLOR="#BC8F8F">&quot;center_y&quot;</FONT></B>, online_center_y);
      rb_eval.set_parameters(online_mu);
      rb_eval.print_parameters();
      rb_eval.rb_solve( rb_eval.get_n_basis_functions() );
  
      <B><FONT COLOR="#A020F0">if</FONT></B>(store_basis_functions)
      {
        eim_rb_eval.read_in_basis_functions(eim_construction,<B><FONT COLOR="#BC8F8F">&quot;eim_data&quot;</FONT></B>);
        rb_eval.read_in_basis_functions(rb_construction,<B><FONT COLOR="#BC8F8F">&quot;rb_data&quot;</FONT></B>);
  
        eim_construction.load_rb_solution();
        rb_construction.load_rb_solution();
  #ifdef LIBMESH_HAVE_EXODUS_API
        ExodusII_IO(mesh).write_equation_systems(<B><FONT COLOR="#BC8F8F">&quot;RB_sol.e&quot;</FONT></B>,equation_systems);
  #endif
      }
    }
  
  }
</pre> 
<a name="output"></a> 
<br><br><br> <h1> The console output of the program: </h1> 
<pre>
make[4]: Entering directory `/net/spark/workspace/roystgnr/libmesh/git/devel/examples/reduced_basis/reduced_basis_ex4'
***************************************************************
* Running Example reduced_basis_ex4:
*  mpirun -np 4 example-devel -online_mode 0 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc
***************************************************************
 
*** Warning, This code is untested, experimental, or likely to see future API changes: ../src/reduced_basis/rb_parametrized.C, line 41, compiled Apr 19 2013 at 11:42:51 ***
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

 EquationSystems
  n_systems()=2
   System #0, "EIM"
    Type "RBConstruction"
    Variables="f_EIM" 
    Finite Element Types="LAGRANGE", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" 
    Approximation Orders="FIRST", "THIRD" 
    n_dofs()=676
    n_local_dofs()=186
    n_constrained_dofs()=0
    n_local_constrained_dofs()=0
    n_vectors()=1
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 8.33136
      Average Off-Processor Bandwidth <= 0.488166
      Maximum  On-Processor Bandwidth <= 11
      Maximum Off-Processor Bandwidth <= 8
    DofMap Constraints
      Number of DoF Constraints = 0
      Number of Node Constraints = 0
   System #1, "RB"
    Type "RBConstruction"
    Variables="u" 
    Finite Element Types="LAGRANGE", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" 
    Approximation Orders="FIRST", "THIRD" 
    n_dofs()=676
    n_local_dofs()=186
    n_constrained_dofs()=100
    n_local_constrained_dofs()=25
    n_vectors()=1
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

Initializing training parameters with deterministic training set...
Parameter center_x: log scaling = 0
Parameter center_y: log scaling = 0


RBConstruction parameters:
system name: EIM
constrained_problem: 0
Nmax: 20
Basis training error tolerance: 0.001
Aq operators attached: 0
Fq functions attached: 0
n_outputs: 0
Number of parameters: 2
Parameter center_x: Min = -1, Max = 1, value = 1
Parameter center_y: Min = -1, Max = 1, value = 1
n_training_samples: 25
single-matrix mode? 0
reuse preconditioner? 1
use a relative error bound in greedy? 0
write out data during basis training? 0
quiet mode? 1


RBEIMConstruction parameters:
best fit type: projection

Initializing parametrized functions in training set...
Completed solve for training sample 1 of 25
Completed solve for training sample 2 of 25
Completed solve for training sample 3 of 25
Completed solve for training sample 4 of 25
Completed solve for training sample 5 of 25
Completed solve for training sample 6 of 25
Completed solve for training sample 7 of 25
Completed solve for training sample 8 of 25
Completed solve for training sample 9 of 25
Completed solve for training sample 10 of 25
Completed solve for training sample 11 of 25
Completed solve for training sample 12 of 25
Completed solve for training sample 13 of 25
Completed solve for training sample 14 of 25
Completed solve for training sample 15 of 25
Completed solve for training sample 16 of 25
Completed solve for training sample 17 of 25
Completed solve for training sample 18 of 25
Completed solve for training sample 19 of 25
Completed solve for training sample 20 of 25
Completed solve for training sample 21 of 25
Completed solve for training sample 22 of 25
Completed solve for training sample 23 of 25
Completed solve for training sample 24 of 25
Completed solve for training sample 25 of 25
Parametrized functions in training set initialized


---- Performing Greedy basis enrichment ----

---- Basis dimension: 0 ----
Performing truth solve at parameter:
center_x: 1.000000e+00
center_y: 1.000000e+00

Enriching the RB space
Updating RB matrices

---- Basis dimension: 1 ----
Performing RB solves on training set
Maximum (absolute) error bound is 1.04028

Performing truth solve at parameter:
center_x: 5.000000e-01
center_y: 5.000000e-01

Enriching the RB space
Updating RB matrices

---- Basis dimension: 2 ----
Performing RB solves on training set
Maximum (absolute) error bound is 1.00428

Performing truth solve at parameter:
center_x: -1.000000e+00
center_y: -1.000000e+00

Enriching the RB space
Updating RB matrices

---- Basis dimension: 3 ----
Performing RB solves on training set
Maximum (absolute) error bound is 1.02817

Performing truth solve at parameter:
center_x: -5.000000e-01
center_y: -5.000000e-01

Enriching the RB space
Updating RB matrices

---- Basis dimension: 4 ----
Performing RB solves on training set
Maximum (absolute) error bound is 1.00298

Performing truth solve at parameter:
center_x: 1.000000e+00
center_y: -1.000000e+00

Enriching the RB space
Updating RB matrices

---- Basis dimension: 5 ----
Performing RB solves on training set
Maximum (absolute) error bound is 1.00292

Performing truth solve at parameter:
center_x: -1.000000e+00
center_y: 1.000000e+00

Enriching the RB space
Updating RB matrices

---- Basis dimension: 6 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.86594

Performing truth solve at parameter:
center_x: 5.000000e-01
center_y: -5.000000e-01

Enriching the RB space
Updating RB matrices

---- Basis dimension: 7 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.816903

Performing truth solve at parameter:
center_x: -5.000000e-01
center_y: 5.000000e-01

Enriching the RB space
Updating RB matrices

---- Basis dimension: 8 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.681622

Performing truth solve at parameter:
center_x: -1.000000e+00
center_y: 0.000000e+00

Enriching the RB space
Updating RB matrices

---- Basis dimension: 9 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.658046

Performing truth solve at parameter:
center_x: 1.000000e+00
center_y: 0.000000e+00

Enriching the RB space
Updating RB matrices

---- Basis dimension: 10 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.606155

Performing truth solve at parameter:
center_x: 0.000000e+00
center_y: -1.000000e+00

Enriching the RB space
Updating RB matrices

---- Basis dimension: 11 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.606138

Performing truth solve at parameter:
center_x: 0.000000e+00
center_y: 1.000000e+00

Enriching the RB space
Updating RB matrices

---- Basis dimension: 12 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.193607

Performing truth solve at parameter:
center_x: 0.000000e+00
center_y: 0.000000e+00

Enriching the RB space
Updating RB matrices

---- Basis dimension: 13 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.184868

Performing truth solve at parameter:
center_x: -1.000000e+00
center_y: 5.000000e-01

Enriching the RB space
Updating RB matrices

---- Basis dimension: 14 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.214039

Performing truth solve at parameter:
center_x: -5.000000e-01
center_y: 1.000000e+00

Enriching the RB space
Updating RB matrices

---- Basis dimension: 15 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.175349

Performing truth solve at parameter:
center_x: 1.000000e+00
center_y: 5.000000e-01

Enriching the RB space
Updating RB matrices

---- Basis dimension: 16 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.163258

Performing truth solve at parameter:
center_x: 5.000000e-01
center_y: 1.000000e+00

Enriching the RB space
Updating RB matrices

---- Basis dimension: 17 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.161335

Performing truth solve at parameter:
center_x: 5.000000e-01
center_y: -1.000000e+00

Enriching the RB space
Updating RB matrices

---- Basis dimension: 18 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.159752

Performing truth solve at parameter:
center_x: 1.000000e+00
center_y: -5.000000e-01

Enriching the RB space
Updating RB matrices

---- Basis dimension: 19 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.123948

Performing truth solve at parameter:
center_x: -1.000000e+00
center_y: -5.000000e-01

Enriching the RB space
Updating RB matrices

---- Basis dimension: 20 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.162888

Maximum number of basis functions reached: Nmax = 20.
Perform one more Greedy iteration for error bounds.
Performing truth solve at parameter:
center_x: -5.000000e-01
center_y: -1.000000e+00

Enriching the RB space
Updating RB matrices

---- Basis dimension: 20 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.162888

Extra Greedy iteration finished.
Initializing training parameters with deterministic training set...
Parameter center_x: log scaling = 0
Parameter center_y: log scaling = 0


RBConstruction parameters:
system name: RB
constrained_problem: 0
Nmax: 15
Basis training error tolerance: 0.001
Aq operators attached: 1
Fq functions attached: 20
n_outputs: 0
Number of parameters: 2
Parameter center_x: Min = -1, Max = 1, value = 1
Parameter center_y: Min = -1, Max = 1, value = 1
n_training_samples: 100
single-matrix mode? 0
reuse preconditioner? 1
use a relative error bound in greedy? 0
write out data during basis training? 0
quiet mode? 1


---- Performing Greedy basis enrichment ----

---- Basis dimension: 0 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.388534

Performing truth solve at parameter:
center_x: 1.111111e-01
center_y: -1.111111e-01

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 1 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.17691

Performing truth solve at parameter:
center_x: -5.555556e-01
center_y: 5.555556e-01

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 2 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.156396

Performing truth solve at parameter:
center_x: 5.555556e-01
center_y: 5.555556e-01

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 3 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.130261

Performing truth solve at parameter:
center_x: -5.555556e-01
center_y: -5.555556e-01

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 4 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.0995149

Performing truth solve at parameter:
center_x: 5.555556e-01
center_y: -5.555556e-01

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 5 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.0478085

Performing truth solve at parameter:
center_x: -1.111111e-01
center_y: -1.000000e+00

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 6 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.0432272

Performing truth solve at parameter:
center_x: 1.111111e-01
center_y: 1.000000e+00

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 7 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.0306345

Performing truth solve at parameter:
center_x: 1.000000e+00
center_y: -1.111111e-01

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 8 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.0268713

Performing truth solve at parameter:
center_x: -1.000000e+00
center_y: 1.111111e-01

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 9 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.0243775

Performing truth solve at parameter:
center_x: -1.000000e+00
center_y: 1.000000e+00

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 10 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.0231995

Performing truth solve at parameter:
center_x: -1.000000e+00
center_y: -1.000000e+00

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 11 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.0207146

Performing truth solve at parameter:
center_x: 1.000000e+00
center_y: 1.000000e+00

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 12 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.0181294

Performing truth solve at parameter:
center_x: 1.000000e+00
center_y: -1.000000e+00

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 13 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.00827867

Performing truth solve at parameter:
center_x: -3.333333e-01
center_y: 1.000000e+00

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 14 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.00726206

Performing truth solve at parameter:
center_x: 5.555556e-01
center_y: -1.000000e+00

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 15 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.00574847

Maximum number of basis functions reached: Nmax = 15

 -------------------------------------------------------------------------------------------------------------------
| Processor id:   0                                                                                                 |
| Num Processors: 4                                                                                                 |
| Time:           Fri Apr 19 11:55:34 2013                                                                          |
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
| libMesh Performance: Alive time=2.19545, Active time=2.13119                                                      |
 -------------------------------------------------------------------------------------------------------------------
| Event                                 nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                                 w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|-------------------------------------------------------------------------------------------------------------------|
|                                                                                                                   |
|                                                                                                                   |
| DofMap                                                                                                            |
|   add_neighbors_to_send_list()        2         0.0015      0.000735    0.0021      0.001058    0.07     0.10     |
|   build_constraint_matrix()           3432      0.0041      0.000001    0.0041      0.000001    0.19     0.19     |
|   build_sparsity()                    2         0.0012      0.000601    0.0033      0.001635    0.06     0.15     |
|   cnstrn_elem_mat_vec()               3432      0.0052      0.000002    0.0052      0.000002    0.25     0.25     |
|   create_dof_constraints()            2         0.0020      0.000984    0.0049      0.002474    0.09     0.23     |
|   distribute_dofs()                   2         0.0030      0.001508    0.0174      0.008693    0.14     0.82     |
|   dof_indices()                       25849     0.1102      0.000004    0.1102      0.000004    5.17     5.17     |
|   prepare_send_list()                 2         0.0000      0.000007    0.0000      0.000007    0.00     0.00     |
|   reinit()                            2         0.0048      0.002407    0.0048      0.002407    0.23     0.23     |
|                                                                                                                   |
| FE                                                                                                                |
|   compute_shape_functions()           22632     0.0876      0.000004    0.0876      0.000004    4.11     4.11     |
|   init_shape_functions()              1242      0.0058      0.000005    0.0058      0.000005    0.27     0.27     |
|   inverse_map()                       15568     0.0395      0.000003    0.0395      0.000003    1.85     1.85     |
|                                                                                                                   |
| FEMap                                                                                                             |
|   compute_affine_map()                22632     0.0507      0.000002    0.0507      0.000002    2.38     2.38     |
|   compute_face_map()                  1104      0.0046      0.000004    0.0104      0.000009    0.22     0.49     |
|   init_face_shape_functions()         46        0.0001      0.000002    0.0001      0.000002    0.00     0.00     |
|   init_reference_to_physical_map()    1242      0.0046      0.000004    0.0046      0.000004    0.22     0.22     |
|                                                                                                                   |
| Mesh                                                                                                              |
|   find_neighbors()                    1         0.0017      0.001710    0.0020      0.001993    0.08     0.09     |
|   renumber_nodes_and_elem()           2         0.0002      0.000081    0.0002      0.000081    0.01     0.01     |
|                                                                                                                   |
| MeshCommunication                                                                                                 |
|   assign_global_indices()             2         0.0157      0.007860    0.0249      0.012469    0.74     1.17     |
|   compute_hilbert_indices()           2         0.0050      0.002521    0.0050      0.002521    0.24     0.24     |
|   find_global_indices()               2         0.0007      0.000365    0.0153      0.007646    0.03     0.72     |
|   parallel_sort()                     2         0.0005      0.000231    0.0092      0.004598    0.02     0.43     |
|                                                                                                                   |
| MeshTools::Generation                                                                                             |
|   build_cube()                        1         0.0005      0.000536    0.0005      0.000536    0.03     0.03     |
|                                                                                                                   |
| MetisPartitioner                                                                                                  |
|   partition()                         1         0.0050      0.005015    0.0082      0.008246    0.24     0.39     |
|                                                                                                                   |
| Parallel                                                                                                          |
|   allgather()                         22        0.0163      0.000741    0.0164      0.000747    0.76     0.77     |
|   barrier()                           2         0.0000      0.000009    0.0000      0.000009    0.00     0.00     |
|   broadcast()                         560       0.0024      0.000004    0.0024      0.000004    0.11     0.11     |
|   max(bool)                           5         0.0000      0.000004    0.0000      0.000004    0.00     0.00     |
|   max(scalar)                         617       0.0099      0.000016    0.0099      0.000016    0.46     0.46     |
|   max(vector)                         42        0.0003      0.000008    0.0010      0.000025    0.02     0.05     |
|   maxloc(scalar)                      58        0.1202      0.002073    0.1202      0.002073    5.64     5.64     |
|   min(bool)                           202       0.0010      0.000005    0.0010      0.000005    0.05     0.05     |
|   min(scalar)                         166       0.0136      0.000082    0.0136      0.000082    0.64     0.64     |
|   min(vector)                         42        0.0004      0.000011    0.0013      0.000031    0.02     0.06     |
|   probe()                             142       0.0041      0.000029    0.0041      0.000029    0.19     0.19     |
|   receive()                           110       0.0004      0.000003    0.0008      0.000007    0.02     0.04     |
|   send()                              86        0.0002      0.000002    0.0002      0.000002    0.01     0.01     |
|   send_receive()                      90        0.0004      0.000004    0.0012      0.000014    0.02     0.06     |
|   sum()                               56        0.0102      0.000182    0.0103      0.000184    0.48     0.48     |
|                                                                                                                   |
| Parallel::Request                                                                                                 |
|   wait()                              86        0.0001      0.000001    0.0001      0.000001    0.00     0.00     |
|                                                                                                                   |
| Partitioner                                                                                                       |
|   set_node_processor_ids()            1         0.0005      0.000524    0.0007      0.000674    0.02     0.03     |
|   set_parent_processor_ids()          1         0.0002      0.000187    0.0002      0.000187    0.01     0.01     |
|                                                                                                                   |
| PetscLinearSolver                                                                                                 |
|   solve()                             75        0.2182      0.002909    0.2182      0.002909    10.24    10.24    |
|                                                                                                                   |
| PointLocatorTree                                                                                                  |
|   init(no master)                     1         0.0013      0.001298    0.0015      0.001547    0.06     0.07     |
|   operator()                          440       0.0126      0.000029    0.0148      0.000034    0.59     0.69     |
|                                                                                                                   |
| RBConstruction                                                                                                    |
|   add_scaled_matrix_and_vector()      23        0.1996      0.008679    0.3434      0.014932    9.37     16.12    |
|   clear()                             3         0.0006      0.000186    0.0006      0.000186    0.03     0.03     |
|   compute_Fq_representor_innerprods() 2         0.0047      0.002338    0.0500      0.025020    0.22     2.35     |
|   compute_max_error_bound()           37        0.0067      0.000182    1.0228      0.027642    0.32     47.99    |
|   enrich_RB_space()                   15        0.0087      0.000580    0.0087      0.000580    0.41     0.41     |
|   train_reduced_basis()               2         0.0095      0.004750    1.4281      0.714063    0.45     67.01    |
|   truth_assembly()                    15        0.0176      0.001174    0.0183      0.001222    0.83     0.86     |
|   truth_solve()                       15        0.0015      0.000102    0.0591      0.003938    0.07     2.77     |
|   update_RB_system_matrices()         36        0.0278      0.000773    0.0278      0.000773    1.31     1.31     |
|   update_residual_terms()             15        0.0299      0.001996    0.0701      0.004673    1.40     3.29     |
|                                                                                                                   |
| RBEIMConstruction                                                                                                 |
|   compute_best_fit_error()            525       0.5501      0.001048    0.6351      0.001210    25.81    29.80    |
|   enrich_RB_space()                   21        0.0430      0.002050    0.1512      0.007199    2.02     7.09     |
|   truth_solve()                       571       0.1651      0.000289    0.3488      0.000611    7.74     16.37    |
|   update_RB_system_matrices()         21        0.0170      0.000808    0.0429      0.002045    0.80     2.01     |
|                                                                                                                   |
| RBEIMEvaluation                                                                                                   |
|   rb_solve()                          434       0.0150      0.000034    0.0150      0.000034    0.70     0.70     |
|   write_offline_data_to_files()       1         0.0004      0.000374    0.0011      0.001116    0.02     0.05     |
|                                                                                                                   |
| RBEvaluation                                                                                                      |
|   clear()                             3         0.0001      0.000033    0.0001      0.000033    0.00     0.00     |
|   compute_residual_dual_norm()        400       0.2404      0.000601    0.2404      0.000601    11.28    11.28    |
|   rb_solve()                          400       0.0157      0.000039    0.2700      0.000675    0.74     12.67    |
|   resize_data_structures()            2         0.0000      0.000024    0.0000      0.000024    0.00     0.00     |
|   write_offline_data_to_files()       2         0.0023      0.001133    0.0023      0.001133    0.11     0.11     |
|   write_out_basis_functions()         2         0.0003      0.000133    0.0402      0.020101    0.01     1.89     |
|   write_out_vectors()                 2         0.0086      0.004292    0.0399      0.019967    0.40     1.87     |
 -------------------------------------------------------------------------------------------------------------------
| Totals:                               102555    2.1312                                          100.00            |
 -------------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example reduced_basis_ex4:
*  mpirun -np 4 example-devel -online_mode 0 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc
***************************************************************
***************************************************************
* Running Example reduced_basis_ex4:
*  mpirun -np 4 example-devel -online_mode 1 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc
***************************************************************
 
*** Warning, This code is untested, experimental, or likely to see future API changes: ../src/reduced_basis/rb_parametrized.C, line 41, compiled Apr 19 2013 at 11:42:51 ***
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

 EquationSystems
  n_systems()=2
   System #0, "EIM"
    Type "RBConstruction"
    Variables="f_EIM" 
    Finite Element Types="LAGRANGE", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" 
    Approximation Orders="FIRST", "THIRD" 
    n_dofs()=676
    n_local_dofs()=186
    n_constrained_dofs()=0
    n_local_constrained_dofs()=0
    n_vectors()=1
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 8.33136
      Average Off-Processor Bandwidth <= 0.488166
      Maximum  On-Processor Bandwidth <= 11
      Maximum Off-Processor Bandwidth <= 8
    DofMap Constraints
      Number of DoF Constraints = 0
      Number of Node Constraints = 0
   System #1, "RB"
    Type "RBConstruction"
    Variables="u" 
    Finite Element Types="LAGRANGE", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" 
    Approximation Orders="FIRST", "THIRD" 
    n_dofs()=676
    n_local_dofs()=186
    n_constrained_dofs()=100
    n_local_constrained_dofs()=25
    n_vectors()=1
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

center_x: -6.000000e-01
center_y: 7.000000e-01


 -------------------------------------------------------------------------------------------------------------------
| Processor id:   0                                                                                                 |
| Num Processors: 4                                                                                                 |
| Time:           Fri Apr 19 11:55:34 2013                                                                          |
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
| libMesh Performance: Alive time=0.21897, Active time=0.196295                                                |
 --------------------------------------------------------------------------------------------------------------
| Event                            nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                            w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|--------------------------------------------------------------------------------------------------------------|
|                                                                                                              |
|                                                                                                              |
| DofMap                                                                                                       |
|   add_neighbors_to_send_list()   2         0.0015      0.000745    0.0021      0.001067    0.76     1.09     |
|   build_sparsity()               2         0.0013      0.000638    0.0034      0.001701    0.65     1.73     |
|   create_dof_constraints()       2         0.0020      0.001010    0.0049      0.002435    1.03     2.48     |
|   distribute_dofs()              2         0.0031      0.001548    0.0095      0.004754    1.58     4.84     |
|   dof_indices()                  1387      0.0063      0.000005    0.0063      0.000005    3.23     3.23     |
|   prepare_send_list()            2         0.0000      0.000007    0.0000      0.000007    0.01     0.01     |
|   reinit()                       2         0.0051      0.002543    0.0051      0.002543    2.59     2.59     |
|                                                                                                              |
| EquationSystems                                                                                              |
|   build_solution_vector()        1         0.0008      0.000816    0.0028      0.002849    0.42     1.45     |
|                                                                                                              |
| ExodusII_IO                                                                                                  |
|   write_nodal_data()             1         0.0935      0.093465    0.0935      0.093465    47.61    47.61    |
|                                                                                                              |
| Mesh                                                                                                         |
|   find_neighbors()               1         0.0018      0.001761    0.0020      0.001970    0.90     1.00     |
|   renumber_nodes_and_elem()      2         0.0002      0.000076    0.0002      0.000076    0.08     0.08     |
|                                                                                                              |
| MeshCommunication                                                                                            |
|   assign_global_indices()        2         0.0270      0.013520    0.0281      0.014054    13.78    14.32    |
|   compute_hilbert_indices()      2         0.0055      0.002743    0.0055      0.002743    2.80     2.80     |
|   find_global_indices()          2         0.0008      0.000384    0.0073      0.003626    0.39     3.69     |
|   parallel_sort()                2         0.0004      0.000217    0.0007      0.000329    0.22     0.34     |
|                                                                                                              |
| MeshOutput                                                                                                   |
|   write_equation_systems()       1         0.0001      0.000068    0.0965      0.096462    0.03     49.14    |
|                                                                                                              |
| MeshTools::Generation                                                                                        |
|   build_cube()                   1         0.0005      0.000530    0.0005      0.000530    0.27     0.27     |
|                                                                                                              |
| MetisPartitioner                                                                                             |
|   partition()                    1         0.0052      0.005158    0.0087      0.008722    2.63     4.44     |
|                                                                                                              |
| Parallel                                                                                                     |
|   allgather()                    22        0.0008      0.000038    0.0009      0.000040    0.42     0.45     |
|   barrier()                      2         0.0069      0.003438    0.0069      0.003438    3.50     3.50     |
|   broadcast()                    26        0.0001      0.000004    0.0001      0.000003    0.05     0.04     |
|   max(bool)                      2         0.0000      0.000003    0.0000      0.000003    0.00     0.00     |
|   max(scalar)                    177       0.0015      0.000009    0.0015      0.000009    0.77     0.77     |
|   max(vector)                    39        0.0004      0.000010    0.0014      0.000035    0.21     0.69     |
|   min(bool)                      201       0.0015      0.000008    0.0015      0.000008    0.78     0.78     |
|   min(scalar)                    166       0.0036      0.000022    0.0036      0.000022    1.83     1.83     |
|   min(vector)                    39        0.0005      0.000014    0.0015      0.000039    0.27     0.77     |
|   probe()                        110       0.0004      0.000004    0.0004      0.000004    0.20     0.20     |
|   receive()                      98        0.0004      0.000004    0.0007      0.000008    0.18     0.37     |
|   send()                         98        0.0003      0.000003    0.0003      0.000003    0.15     0.15     |
|   send_receive()                 90        0.0004      0.000005    0.0013      0.000015    0.21     0.69     |
|   sum()                          48        0.0006      0.000011    0.0009      0.000018    0.28     0.43     |
|                                                                                                              |
| Parallel::Request                                                                                            |
|   wait()                         98        0.0001      0.000001    0.0001      0.000001    0.07     0.07     |
|                                                                                                              |
| Partitioner                                                                                                  |
|   set_node_processor_ids()       1         0.0005      0.000520    0.0012      0.001178    0.26     0.60     |
|   set_parent_processor_ids()     1         0.0002      0.000164    0.0002      0.000164    0.08     0.08     |
|                                                                                                              |
| RBConstruction                                                                                               |
|   clear()                        3         0.0004      0.000135    0.0004      0.000135    0.21     0.21     |
|   load_rb_solution()             2         0.0002      0.000118    0.0002      0.000118    0.12     0.12     |
|                                                                                                              |
| RBEIMEvaluation                                                                                              |
|   rb_solve()                     1         0.0088      0.008784    0.0088      0.008784    4.47     4.47     |
|   read_offline_data_from_files() 1         0.0002      0.000213    0.0006      0.000636    0.11     0.32     |
|                                                                                                              |
| RBEvaluation                                                                                                 |
|   clear()                        3         0.0001      0.000044    0.0001      0.000044    0.07     0.07     |
|   compute_residual_dual_norm()   1         0.0010      0.000968    0.0010      0.000968    0.49     0.49     |
|   rb_solve()                     1         0.0002      0.000174    0.0099      0.009926    0.09     5.06     |
|   read_in_basis_functions()      2         0.0000      0.000013    0.0467      0.023363    0.01     23.80    |
|   read_in_vectors()              2         0.0111      0.005526    0.0467      0.023349    5.63     23.79    |
|   read_offline_data_from_files() 2         0.0010      0.000509    0.0011      0.000549    0.52     0.56     |
|   resize_data_structures()       2         0.0001      0.000040    0.0001      0.000040    0.04     0.04     |
 --------------------------------------------------------------------------------------------------------------
| Totals:                          2653      0.1963                                          100.00            |
 --------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example reduced_basis_ex4:
*  mpirun -np 4 example-devel -online_mode 1 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc
***************************************************************
make[4]: Leaving directory `/net/spark/workspace/roystgnr/libmesh/git/devel/examples/reduced_basis/reduced_basis_ex4'
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
