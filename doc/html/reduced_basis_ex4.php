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
        
          SimpleEIMEvaluation()
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
                                 const std::string& name,
                                 const unsigned int number)
          : Parent(es, name, number)
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
          SimpleRBEvaluation()
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
<h1>Reduced Basis Example 4 - Empirical Interpolation Method</h1>


<br><br>In this example problem we develop a reduced basis approximation for a parametrized
PDE that has "non-affine" parameter dependence. This requires the use of the
Empirical Interpolation Method (EIM).

<br><br>We first use EIM to construct an affine approximation to the non-affine term,
which is a parametrized function that is a Gaussian with "center" defined
by the two parameters (mu_1,mu_2) \in [-1,1]^2. We then employ this EIM
approximation in order to generate a reduced basis approximation for the
parametrized PDE: -0.05 * Laplacian(u) = f(mu_1,mu_2), with zero Dirichlet
boundary conditions.


<br><br>Basic include file needed for the mesh functionality.
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
Create a mesh (just a simple square)
</div>

<div class ="fragment">
<pre>
          Mesh mesh (dim);
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
  
    SimpleEIMEvaluation()
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
                           <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp; name,
                           <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> number)
    : Parent(es, name, number)
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
    SimpleRBEvaluation()
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
  
    Mesh mesh (dim);
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
  
    SimpleRBEvaluation rb_eval;
  
    SimpleEIMEvaluation eim_rb_eval;
    
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
***************************************************************
* Running Example reduced_basis_ex4:
*  mpirun -np 12 example-devel -online_mode 0 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
***************************************************************
 
*** Warning, This code is untested, experimental, or likely to see future API changes: src/reduced_basis/rb_parametrized.C, line 41, compiled Jan 31 2013 at 21:51:32 ***
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

 EquationSystems
  n_systems()=2
   System #0, "EIM"
    Type "RBConstruction"
    Variables="f_EIM" 
    Finite Element Types="LAGRANGE", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" 
    Approximation Orders="FIRST", "THIRD" 
    n_dofs()=676
    n_local_dofs()=68
    n_constrained_dofs()=0
    n_local_constrained_dofs()=0
    n_vectors()=1
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 7.99852
      Average Off-Processor Bandwidth <= 1.20118
      Maximum  On-Processor Bandwidth <= 11
      Maximum Off-Processor Bandwidth <= 9
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
    n_local_dofs()=68
    n_constrained_dofs()=100
    n_local_constrained_dofs()=16
    n_vectors()=1
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
center_x: -5.000000e-01
center_y: 5.000000e-01

Enriching the RB space
Updating RB matrices

---- Basis dimension: 7 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.816903

Performing truth solve at parameter:
center_x: 5.000000e-01
center_y: -5.000000e-01

Enriching the RB space
Updating RB matrices

---- Basis dimension: 8 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.681622

Performing truth solve at parameter:
center_x: 1.000000e+00
center_y: 0.000000e+00

Enriching the RB space
Updating RB matrices

---- Basis dimension: 9 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.658046

Performing truth solve at parameter:
center_x: -1.000000e+00
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
center_x: 5.000000e-01
center_y: -1.000000e+00

Enriching the RB space
Updating RB matrices

---- Basis dimension: 14 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.214039

Performing truth solve at parameter:
center_x: 1.000000e+00
center_y: -5.000000e-01

Enriching the RB space
Updating RB matrices

---- Basis dimension: 15 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.175349

Performing truth solve at parameter:
center_x: -1.000000e+00
center_y: -5.000000e-01

Enriching the RB space
Updating RB matrices

---- Basis dimension: 16 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.163258

Performing truth solve at parameter:
center_x: -5.000000e-01
center_y: -1.000000e+00

Enriching the RB space
Updating RB matrices

---- Basis dimension: 17 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.161335

Performing truth solve at parameter:
center_x: 5.000000e-01
center_y: 1.000000e+00

Enriching the RB space
Updating RB matrices

---- Basis dimension: 18 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.159752

Performing truth solve at parameter:
center_x: 1.000000e+00
center_y: 5.000000e-01

Enriching the RB space
Updating RB matrices

---- Basis dimension: 19 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.123948

Performing truth solve at parameter:
center_x: -5.000000e-01
center_y: 1.000000e+00

Enriching the RB space
Updating RB matrices

---- Basis dimension: 20 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.162888

Maximum number of basis functions reached: Nmax = 20.
Perform one more Greedy iteration for error bounds.
Performing truth solve at parameter:
center_x: -1.000000e+00
center_y: 5.000000e-01

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
Maximum (absolute) error bound is 0.388376

Performing truth solve at parameter:
center_x: -1.111111e-01
center_y: 1.111111e-01

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 1 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.176776

Performing truth solve at parameter:
center_x: 5.555556e-01
center_y: -5.555556e-01

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 2 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.15644

Performing truth solve at parameter:
center_x: 5.555556e-01
center_y: 5.555556e-01

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 3 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.130994

Performing truth solve at parameter:
center_x: -5.555556e-01
center_y: -5.555556e-01

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 4 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.0984295

Performing truth solve at parameter:
center_x: -5.555556e-01
center_y: 5.555556e-01

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 5 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.0475053

Performing truth solve at parameter:
center_x: -1.000000e+00
center_y: 1.111111e-01

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 6 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.0432771

Performing truth solve at parameter:
center_x: 1.000000e+00
center_y: -1.111111e-01

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 7 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.0310976

Performing truth solve at parameter:
center_x: 1.111111e-01
center_y: 1.000000e+00

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 8 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.0264776

Performing truth solve at parameter:
center_x: -1.111111e-01
center_y: -1.000000e+00

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 9 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.0242248

Performing truth solve at parameter:
center_x: 1.000000e+00
center_y: -1.000000e+00

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 10 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.0228554

Performing truth solve at parameter:
center_x: -1.000000e+00
center_y: -1.000000e+00

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 11 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.0206759

Performing truth solve at parameter:
center_x: 1.000000e+00
center_y: 1.000000e+00

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 12 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.0188747

Performing truth solve at parameter:
center_x: -1.000000e+00
center_y: 1.000000e+00

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 13 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.00819468

Performing truth solve at parameter:
center_x: 3.333333e-01
center_y: -1.000000e+00

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 14 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.00761207

Performing truth solve at parameter:
center_x: -3.333333e-01
center_y: 1.000000e+00

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 15 ----
Performing RB solves on training set
Maximum (absolute) error bound is 0.00555548

Maximum number of basis functions reached: Nmax = 15
************************************************************************************************************************
***             WIDEN YOUR WINDOW TO 120 CHARACTERS.  Use 'enscript -r -fCourier9' to print this document            ***
************************************************************************************************************************

---------------------------------------------- PETSc Performance Summary: ----------------------------------------------

/workspace/libmesh/examples/reduced_basis/reduced_basis_ex4/.libs/lt-example-devel on a intel-12. named hbar.ices.utexas.edu with 12 processors, by benkirk Thu Jan 31 22:18:08 2013
Using Petsc Release Version 3.3.0, Patch 2, Fri Jul 13 15:42:00 CDT 2012 

                         Max       Max/Min        Avg      Total 
Time (sec):           3.712e+00      1.00000   3.712e+00
Objects:              9.080e+02      1.00000   9.080e+02
Flops:                3.442e+07      1.81602   2.798e+07  3.357e+08
Flops/sec:            9.274e+06      1.81602   7.537e+06  9.044e+07
MPI Messages:         6.874e+04      3.02107   4.567e+04  5.480e+05
MPI Message Lengths:  3.778e+06      2.28532   5.603e+01  3.070e+07
MPI Reductions:       3.640e+04      1.00000

Flop counting convention: 1 flop = 1 real number operation of type (multiply/divide/add/subtract)
                            e.g., VecAXPY() for real vectors of length N --> 2N flops
                            and VecAXPY() for complex vectors of length N --> 8N flops

Summary of Stages:   ----- Time ------  ----- Flops -----  --- Messages ---  -- Message Lengths --  -- Reductions --
                        Avg     %Total     Avg     %Total   counts   %Total     Avg         %Total   counts   %Total 
 0:      Main Stage: 3.7118e+00 100.0%  3.3571e+08 100.0%  5.480e+05 100.0%  5.603e+01      100.0%  3.640e+04 100.0% 

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

VecDot              7375 1.0 3.4578e-01 2.2 9.96e+05 1.8 0.0e+00 0.0e+00 7.4e+03  7  3  0  0 20   7  3  0  0 20    29
VecMDot             3271 1.0 3.1318e-02 1.2 6.13e+06 1.8 0.0e+00 0.0e+00 3.3e+03  1 18  0  0  9   1 18  0  0  9  1942
VecNorm             3996 1.0 6.1615e-02 1.2 4.72e+05 1.8 0.0e+00 0.0e+00 4.0e+03  2  1  0  0 11   2  1  0  0 11    76
VecScale            3747 1.0 2.0204e-03 1.3 2.55e+05 1.8 0.0e+00 0.0e+00 0.0e+00  0  1  0  0  0   0  1  0  0  0  1254
VecCopy             2084 1.0 1.2670e-03 1.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecSet             11325 1.0 4.8068e-03 1.3 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecAXPY             6635 1.0 1.3803e-02 1.1 9.00e+05 1.8 0.0e+00 0.0e+00 0.0e+00  0  3  0  0  0   0  3  0  0  0   648
VecMAXPY            3396 1.0 4.9820e-03 1.6 6.62e+06 1.8 0.0e+00 0.0e+00 0.0e+00  0 20  0  0  0   0 20  0  0  0 13205
VecAssemblyBegin    2159 1.0 2.8868e-01 1.6 0.00e+00 0.0 2.2e+03 8.1e+01 6.5e+03  6  0  0  1 18   6  0  0  1 18     0
VecAssemblyEnd      2159 1.0 1.4429e-03 1.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecScatterBegin    11209 1.0 4.1280e-02 1.7 0.00e+00 0.0 5.4e+05 5.6e+01 0.0e+00  1  0 99 99  0   1  0 99 99  0     0
VecScatterEnd      11209 1.0 2.3566e-01 8.9 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  4  0  0  0  0   4  0  0  0  0     0
VecNormalize        3396 1.0 2.4838e-02 1.1 6.93e+05 1.8 0.0e+00 0.0e+00 3.4e+03  1  2  0  0  9   1  2  0  0  9   277
MatMult             3396 1.0 3.0067e-02 1.4 3.59e+06 1.6 1.6e+05 5.4e+01 0.0e+00  1 11 30 29  0   1 11 30 29  0  1228
MatMultAdd          6900 1.0 2.2571e-01 4.6 7.76e+06 1.6 3.3e+05 5.4e+01 0.0e+00  4 24 60 58  0   4 24 60 58  0   353
MatSolve            3471 1.0 1.0747e-02 2.3 7.48e+06 2.2 0.0e+00 0.0e+00 0.0e+00  0 19  0  0  0   0 19  0  0  0  6005
MatLUFactorNum        32 1.0 1.3368e-03 2.1 2.31e+05 2.5 0.0e+00 0.0e+00 0.0e+00  0  1  0  0  0   0  1  0  0  0  1407
MatILUFactorSym       32 1.0 2.6917e-03 1.8 0.00e+00 0.0 0.0e+00 0.0e+00 9.6e+01  0  0  0  0  0   0  0  0  0  0     0
MatAssemblyBegin    7119 1.0 4.5067e-01 1.1 0.00e+00 0.0 2.2e+02 2.5e+02 1.4e+04 11  0  0  0 39  11  0  0  0 39     0
MatAssemblyEnd      7119 1.0 2.5899e-02 1.2 0.00e+00 0.0 3.5e+03 1.6e+01 3.0e+02  1  0  1  0  1   1  0  1  0  1     0
MatGetRow           4488 1.8 8.3470e-04 1.8 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatGetRowIJ           32 1.0 1.8358e-05 3.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatGetOrdering        32 1.0 8.9741e-04 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 6.4e+01  0  0  0  0  0   0  0  0  0  0     0
MatZeroEntries        46 1.0 1.0800e-04 1.9 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatAXPY               33 1.0 1.8420e-02 1.0 0.00e+00 0.0 3.2e+03 1.6e+01 5.3e+02  0  0  1  0  1   0  0  1  0  1     0
KSPGMRESOrthog      3271 1.0 3.8336e-02 1.1 1.23e+07 1.8 0.0e+00 0.0e+00 3.3e+03  1 36  0  0  9   1 36  0  0  9  3187
KSPSetUp             107 1.0 5.9795e-04 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPSolve              75 1.0 1.5163e-01 1.0 2.48e+07 1.9 1.6e+05 5.4e+01 6.9e+03  4 71 30 29 19   4 71 30 29 19  1565
PCSetUp               64 1.0 1.1992e-02 1.2 2.31e+05 2.5 0.0e+00 0.0e+00 1.6e+02  0  1  0  0  0   0  1  0  0  0   157
PCSetUpOnBlocks       75 1.0 1.0059e-02 1.2 2.31e+05 2.5 0.0e+00 0.0e+00 1.6e+02  0  1  0  0  0   0  1  0  0  0   187
PCApply             3471 1.0 4.0677e-02 1.2 7.48e+06 2.2 0.0e+00 0.0e+00 0.0e+00  1 19  0  0  0   1 19  0  0  0  1586
------------------------------------------------------------------------------------------------------------------------

Memory usage is given in bytes:

Object Type          Creations   Destructions     Memory  Descendants' Mem.
Reports information only for process 0.

--- Event Stage 0: Main Stage

              Vector   451            451       867528     0
      Vector Scatter    61             61        63196     0
           Index Set   218            218       175128     0
   IS L to G Mapping    23             23        12972     0
              Matrix   146            146      1141880     0
       Krylov Solver     4              4        38720     0
      Preconditioner     4              4         3568     0
              Viewer     1              0            0     0
========================================================================================================================
Average time to get PetscTime(): 0
Average time for MPI_Barrier(): 4.62532e-06
Average time for zero size MPI_Send(): 1.32521e-05
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
| Time:           Thu Jan 31 22:18:08 2013                                                                             |
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
| libMesh Performance: Alive time=3.89743, Active time=3.65681                                                      |
 -------------------------------------------------------------------------------------------------------------------
| Event                                 nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                                 w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|-------------------------------------------------------------------------------------------------------------------|
|                                                                                                                   |
|                                                                                                                   |
| DofMap                                                                                                            |
|   add_neighbors_to_send_list()        2         0.0062      0.003114    0.0100      0.005013    0.17     0.27     |
|   build_constraint_matrix()           1122      0.0115      0.000010    0.0115      0.000010    0.31     0.31     |
|   build_sparsity()                    2         0.0043      0.002171    0.0131      0.006570    0.12     0.36     |
|   cnstrn_elem_mat_vec()               1122      0.0097      0.000009    0.0097      0.000009    0.27     0.27     |
|   create_dof_constraints()            2         0.0128      0.006393    0.0474      0.023677    0.35     1.29     |
|   distribute_dofs()                   2         0.0223      0.011149    0.0836      0.041783    0.61     2.29     |
|   dof_indices()                       8897      0.5055      0.000057    0.5055      0.000057    13.82    13.82    |
|   prepare_send_list()                 2         0.0001      0.000035    0.0001      0.000035    0.00     0.00     |
|   reinit()                            2         0.0576      0.028809    0.0576      0.028809    1.58     1.58     |
|                                                                                                                   |
| FE                                                                                                                |
|   compute_shape_functions()           7728      0.1885      0.000024    0.1885      0.000024    5.16     5.16     |
|   init_shape_functions()              828       0.0308      0.000037    0.0308      0.000037    0.84     0.84     |
|   inverse_map()                       6340      0.0649      0.000010    0.0649      0.000010    1.77     1.77     |
|                                                                                                                   |
| FEMap                                                                                                             |
|   compute_affine_map()                7728      0.0794      0.000010    0.0794      0.000010    2.17     2.17     |
|   compute_face_map()                  690       0.0148      0.000021    0.0299      0.000043    0.40     0.82     |
|   init_face_shape_functions()         46        0.0008      0.000017    0.0008      0.000017    0.02     0.02     |
|   init_reference_to_physical_map()    828       0.0182      0.000022    0.0182      0.000022    0.50     0.50     |
|                                                                                                                   |
| Mesh                                                                                                              |
|   find_neighbors()                    1         0.0137      0.013692    0.0148      0.014839    0.37     0.41     |
|   renumber_nodes_and_elem()           2         0.0006      0.000287    0.0006      0.000287    0.02     0.02     |
|                                                                                                                   |
| MeshCommunication                                                                                                 |
|   assign_global_indices()             2         0.0555      0.027760    0.0622      0.031116    1.52     1.70     |
|   compute_hilbert_indices()           2         0.0107      0.005345    0.0107      0.005345    0.29     0.29     |
|   find_global_indices()               2         0.0044      0.002187    0.0200      0.010007    0.12     0.55     |
|   parallel_sort()                     2         0.0031      0.001534    0.0037      0.001827    0.08     0.10     |
|                                                                                                                   |
| MeshTools::Generation                                                                                             |
|   build_cube()                        1         0.0033      0.003255    0.0033      0.003255    0.09     0.09     |
|                                                                                                                   |
| MetisPartitioner                                                                                                  |
|   partition()                         1         0.0492      0.049191    0.0585      0.058473    1.35     1.60     |
|                                                                                                                   |
| Parallel                                                                                                          |
|   allgather()                         22        0.0018      0.000081    0.0019      0.000085    0.05     0.05     |
|   barrier()                           2         0.0001      0.000034    0.0001      0.000034    0.00     0.00     |
|   broadcast()                         560       0.0035      0.000006    0.0035      0.000006    0.09     0.09     |
|   max(bool)                           5         0.0000      0.000006    0.0000      0.000006    0.00     0.00     |
|   max(scalar)                         617       0.0625      0.000101    0.0625      0.000101    1.71     1.71     |
|   max(vector)                         42        0.0006      0.000014    0.0014      0.000034    0.02     0.04     |
|   maxloc(scalar)                      58        0.1044      0.001799    0.1044      0.001799    2.85     2.85     |
|   min(bool)                           202       0.0013      0.000006    0.0013      0.000006    0.04     0.04     |
|   min(scalar)                         166       0.0115      0.000069    0.0115      0.000069    0.31     0.31     |
|   min(vector)                         42        0.0007      0.000017    0.0016      0.000039    0.02     0.04     |
|   probe()                             478       0.0029      0.000006    0.0029      0.000006    0.08     0.08     |
|   receive()                           382       0.0044      0.000011    0.0069      0.000018    0.12     0.19     |
|   send()                              294       0.0009      0.000003    0.0009      0.000003    0.03     0.03     |
|   send_receive()                      298       0.0023      0.000008    0.0099      0.000033    0.06     0.27     |
|   sum()                               56        0.0010      0.000018    0.0024      0.000043    0.03     0.07     |
|                                                                                                                   |
| Parallel::Request                                                                                                 |
|   wait()                              294       0.0006      0.000002    0.0006      0.000002    0.02     0.02     |
|                                                                                                                   |
| Partitioner                                                                                                       |
|   set_node_processor_ids()            1         0.0016      0.001567    0.0027      0.002725    0.04     0.07     |
|   set_parent_processor_ids()          1         0.0011      0.001143    0.0011      0.001143    0.03     0.03     |
|                                                                                                                   |
| PetscLinearSolver                                                                                                 |
|   solve()                             75        0.1795      0.002393    0.1795      0.002393    4.91     4.91     |
|                                                                                                                   |
| PointLocatorTree                                                                                                  |
|   init(no master)                     1         0.0065      0.006497    0.0069      0.006877    0.18     0.19     |
|   operator()                          440       0.0346      0.000079    0.0400      0.000091    0.95     1.09     |
|                                                                                                                   |
| RBConstruction                                                                                                    |
|   add_scaled_matrix_and_vector()      23        0.1530      0.006653    0.5971      0.025960    4.18     16.33    |
|   clear()                             3         0.0014      0.000478    0.0014      0.000478    0.04     0.04     |
|   compute_Fq_representor_innerprods() 2         0.0059      0.002941    0.0523      0.026160    0.16     1.43     |
|   compute_max_error_bound()           37        0.0143      0.000388    1.5333      0.041440    0.39     41.93    |
|   enrich_RB_space()                   15        0.0105      0.000702    0.0105      0.000702    0.29     0.29     |
|   train_reduced_basis()               2         0.0116      0.005799    2.3018      1.150877    0.32     62.94    |
|   truth_assembly()                    15        0.0346      0.002308    0.0364      0.002425    0.95     0.99     |
|   truth_solve()                       15        0.0018      0.000118    0.0779      0.005196    0.05     2.13     |
|   update_RB_system_matrices()         36        0.0428      0.001190    0.0428      0.001190    1.17     1.17     |
|   update_residual_terms()             15        0.0319      0.002126    0.0718      0.004785    0.87     1.96     |
|                                                                                                                   |
| RBEIMConstruction                                                                                                 |
|   compute_best_fit_error()            525       0.8898      0.001695    1.0499      0.002000    24.33    28.71    |
|   enrich_RB_space()                   21        0.0564      0.002683    0.4241      0.020194    1.54     11.60    |
|   truth_solve()                       571       0.2464      0.000432    0.5349      0.000937    6.74     14.63    |
|   update_RB_system_matrices()         21        0.0316      0.001506    0.1020      0.004856    0.87     2.79     |
|                                                                                                                   |
| RBEIMEvaluation                                                                                                   |
|   rb_solve()                          178       0.0164      0.000092    0.0164      0.000092    0.45     0.45     |
|   write_offline_data_to_files()       1         0.0002      0.000232    0.0010      0.000987    0.01     0.03     |
|                                                                                                                   |
| RBEvaluation                                                                                                      |
|   clear()                             3         0.0002      0.000054    0.0002      0.000054    0.00     0.00     |
|   compute_residual_dual_norm()        144       0.4209      0.002923    0.4209      0.002923    11.51    11.51    |
|   rb_solve()                          144       0.0257      0.000179    0.4604      0.003197    0.70     12.59    |
|   resize_data_structures()            2         0.0005      0.000242    0.0005      0.000242    0.01     0.01     |
|   write_offline_data_to_files()       2         0.0020      0.000992    0.0020      0.000992    0.05     0.05     |
|   write_out_basis_functions()         2         0.0001      0.000044    0.1336      0.066800    0.00     3.65     |
|   write_out_vectors()                 2         0.0694      0.034687    0.1335      0.066756    1.90     3.65     |
 -------------------------------------------------------------------------------------------------------------------
| Totals:                               41169     3.6568                                          100.00            |
 -------------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example reduced_basis_ex4:
*  mpirun -np 12 example-devel -online_mode 0 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
***************************************************************
***************************************************************
* Running Example reduced_basis_ex4:
*  mpirun -np 12 example-devel -online_mode 1 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
***************************************************************
 
*** Warning, This code is untested, experimental, or likely to see future API changes: src/reduced_basis/rb_parametrized.C, line 41, compiled Jan 31 2013 at 21:51:32 ***
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

 EquationSystems
  n_systems()=2
   System #0, "EIM"
    Type "RBConstruction"
    Variables="f_EIM" 
    Finite Element Types="LAGRANGE", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" 
    Approximation Orders="FIRST", "THIRD" 
    n_dofs()=676
    n_local_dofs()=68
    n_constrained_dofs()=0
    n_local_constrained_dofs()=0
    n_vectors()=1
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 7.99852
      Average Off-Processor Bandwidth <= 1.20118
      Maximum  On-Processor Bandwidth <= 11
      Maximum Off-Processor Bandwidth <= 9
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
    n_local_dofs()=68
    n_constrained_dofs()=100
    n_local_constrained_dofs()=16
    n_vectors()=1
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

center_x: -6.000000e-01
center_y: 7.000000e-01

************************************************************************************************************************
***             WIDEN YOUR WINDOW TO 120 CHARACTERS.  Use 'enscript -r -fCourier9' to print this document            ***
************************************************************************************************************************

---------------------------------------------- PETSc Performance Summary: ----------------------------------------------

/workspace/libmesh/examples/reduced_basis/reduced_basis_ex4/.libs/lt-example-devel on a intel-12. named hbar.ices.utexas.edu with 12 processors, by benkirk Thu Jan 31 22:18:09 2013
Using Petsc Release Version 3.3.0, Patch 2, Fri Jul 13 15:42:00 CDT 2012 

                         Max       Max/Min        Avg      Total 
Time (sec):           6.237e-01      1.00855   6.219e-01
Objects:              6.000e+01      1.00000   6.000e+01
Flops:                4.760e+03      1.78947   3.943e+03  4.732e+04
Flops/sec:            7.632e+03      1.77484   6.340e+03  7.608e+04
MPI Messages:         3.900e+01      3.25000   2.500e+01  3.000e+02
MPI Message Lengths:  1.548e+03      2.37423   4.165e+01  1.250e+04
MPI Reductions:       1.810e+02      1.00000

Flop counting convention: 1 flop = 1 real number operation of type (multiply/divide/add/subtract)
                            e.g., VecAXPY() for real vectors of length N --> 2N flops
                            and VecAXPY() for complex vectors of length N --> 8N flops

Summary of Stages:   ----- Time ------  ----- Flops -----  --- Messages ---  -- Message Lengths --  -- Reductions --
                        Avg     %Total     Avg     %Total   counts   %Total     Avg         %Total   counts   %Total 
 0:      Main Stage: 6.2190e-01 100.0%  4.7320e+04 100.0%  3.000e+02 100.0%  4.165e+01      100.0%  1.800e+02  99.4% 

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

VecCopy                2 1.0 1.1921e-05 6.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecSet                45 1.0 2.6941e-05 1.7 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecAXPY               35 1.0 5.3103e-03105.1 4.76e+03 1.8 0.0e+00 0.0e+00 0.0e+00  0100  0  0  0   0100  0  0  0     9
VecAssemblyBegin      37 1.0 1.5505e-02 1.5 0.00e+00 0.0 0.0e+00 0.0e+00 1.1e+02  2  0  0  0 61   2  0  0  0 62     0
VecAssemblyEnd        37 1.0 4.1723e-05 1.7 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecScatterBegin        2 1.0 2.1291e-04 4.0 0.00e+00 0.0 1.0e+02 8.1e+01 0.0e+00  0  0 33 65  0   0  0 33 65  0     0
VecScatterEnd          2 1.0 1.9193e-0426.8 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatZeroEntries         4 1.0 1.5974e-05 1.5 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
------------------------------------------------------------------------------------------------------------------------

Memory usage is given in bytes:

Object Type          Creations   Destructions     Memory  Descendants' Mem.
Reports information only for process 0.

--- Event Stage 0: Main Stage

              Vector    45             45        90904     0
      Vector Scatter     2              2         2072     0
           Index Set     4              4         3120     0
   IS L to G Mapping     2              2         1128     0
              Matrix     6              6        33184     0
              Viewer     1              0            0     0
========================================================================================================================
Average time to get PetscTime(): 0
Average time for MPI_Barrier(): 8.39233e-06
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
| Time:           Thu Jan 31 22:18:09 2013                                                                             |
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
| libMesh Performance: Alive time=0.64686, Active time=0.596457                                                |
 --------------------------------------------------------------------------------------------------------------
| Event                            nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                            w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|--------------------------------------------------------------------------------------------------------------|
|                                                                                                              |
|                                                                                                              |
| DofMap                                                                                                       |
|   add_neighbors_to_send_list()   2         0.0063      0.003166    0.0102      0.005088    1.06     1.71     |
|   build_sparsity()               2         0.0044      0.002196    0.0129      0.006464    0.74     2.17     |
|   create_dof_constraints()       2         0.0124      0.006220    0.0470      0.023511    2.09     7.88     |
|   distribute_dofs()              2         0.0220      0.011020    0.1160      0.057981    3.70     19.44    |
|   dof_indices()                  899       0.0492      0.000055    0.0492      0.000055    8.25     8.25     |
|   prepare_send_list()            2         0.0001      0.000032    0.0001      0.000032    0.01     0.01     |
|   reinit()                       2         0.0468      0.023393    0.0468      0.023393    7.84     7.84     |
|                                                                                                              |
| EquationSystems                                                                                              |
|   build_solution_vector()        1         0.0014      0.001429    0.0185      0.018548    0.24     3.11     |
|                                                                                                              |
| ExodusII_IO                                                                                                  |
|   write_nodal_data()             1         0.0070      0.007034    0.0070      0.007034    1.18     1.18     |
|                                                                                                              |
| Mesh                                                                                                         |
|   find_neighbors()               1         0.0137      0.013687    0.0282      0.028157    2.29     4.72     |
|   renumber_nodes_and_elem()      2         0.0005      0.000272    0.0005      0.000272    0.09     0.09     |
|                                                                                                              |
| MeshCommunication                                                                                            |
|   assign_global_indices()        2         0.0671      0.033570    0.0871      0.043543    11.26    14.60    |
|   compute_hilbert_indices()      2         0.0108      0.005404    0.0108      0.005404    1.81     1.81     |
|   find_global_indices()          2         0.0042      0.002119    0.0378      0.018916    0.71     6.34     |
|   parallel_sort()                2         0.0043      0.002127    0.0170      0.008480    0.71     2.84     |
|                                                                                                              |
| MeshOutput                                                                                                   |
|   write_equation_systems()       1         0.0001      0.000096    0.0266      0.026648    0.02     4.47     |
|                                                                                                              |
| MeshTools::Generation                                                                                        |
|   build_cube()                   1         0.0031      0.003069    0.0031      0.003069    0.51     0.51     |
|                                                                                                              |
| MetisPartitioner                                                                                             |
|   partition()                    1         0.0492      0.049177    0.0704      0.070411    8.24     11.80    |
|                                                                                                              |
| Parallel                                                                                                     |
|   allgather()                    22        0.0375      0.001703    0.0379      0.001725    6.28     6.36     |
|   barrier()                      2         0.0019      0.000949    0.0019      0.000949    0.32     0.32     |
|   broadcast()                    26        0.0004      0.000015    0.0003      0.000012    0.06     0.05     |
|   max(bool)                      2         0.0000      0.000007    0.0000      0.000007    0.00     0.00     |
|   max(scalar)                    177       0.0068      0.000038    0.0068      0.000038    1.13     1.13     |
|   max(vector)                    39        0.0015      0.000039    0.0052      0.000134    0.25     0.88     |
|   min(bool)                      201       0.0087      0.000043    0.0087      0.000043    1.45     1.45     |
|   min(scalar)                    166       0.1153      0.000694    0.1153      0.000694    19.32    19.32    |
|   min(vector)                    39        0.0016      0.000041    0.0070      0.000178    0.27     1.17     |
|   probe()                        382       0.0205      0.000054    0.0205      0.000054    3.43     3.43     |
|   receive()                      338       0.0022      0.000007    0.0222      0.000066    0.37     3.73     |
|   send()                         338       0.0011      0.000003    0.0011      0.000003    0.19     0.19     |
|   send_receive()                 298       0.0023      0.000008    0.0257      0.000086    0.38     4.31     |
|   sum()                          48        0.0144      0.000301    0.0305      0.000635    2.42     5.11     |
|                                                                                                              |
| Parallel::Request                                                                                            |
|   wait()                         338       0.0007      0.000002    0.0007      0.000002    0.11     0.11     |
|                                                                                                              |
| Partitioner                                                                                                  |
|   set_node_processor_ids()       1         0.0015      0.001471    0.0529      0.052857    0.25     8.86     |
|   set_parent_processor_ids()     1         0.0012      0.001166    0.0012      0.001166    0.20     0.20     |
|                                                                                                              |
| RBConstruction                                                                                               |
|   clear()                        3         0.0007      0.000219    0.0007      0.000219    0.11     0.11     |
|   load_rb_solution()             2         0.0063      0.003128    0.0063      0.003128    1.05     1.05     |
|                                                                                                              |
| RBEIMEvaluation                                                                                              |
|   rb_solve()                     1         0.0112      0.011224    0.0112      0.011224    1.88     1.88     |
|   read_offline_data_from_files() 1         0.0002      0.000247    0.0011      0.001144    0.04     0.19     |
|                                                                                                              |
| RBEvaluation                                                                                                 |
|   clear()                        3         0.0002      0.000063    0.0002      0.000063    0.03     0.03     |
|   compute_residual_dual_norm()   1         0.0037      0.003678    0.0037      0.003678    0.62     0.62     |
|   rb_solve()                     1         0.0004      0.000385    0.0153      0.015287    0.06     2.56     |
|   read_in_basis_functions()      2         0.0000      0.000024    0.1519      0.075941    0.01     25.46    |
|   read_in_vectors()              2         0.0515      0.025765    0.1518      0.075917    8.64     25.46    |
|   read_offline_data_from_files() 2         0.0016      0.000789    0.0021      0.001066    0.26     0.36     |
|   resize_data_structures()       2         0.0006      0.000276    0.0006      0.000276    0.09     0.09     |
 --------------------------------------------------------------------------------------------------------------
| Totals:                          3365      0.5965                                          100.00            |
 --------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example reduced_basis_ex4:
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
