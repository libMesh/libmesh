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
rbOOmit includes
</div>

<div class ="fragment">
<pre>
        #include "libmesh/rb_parameters.h"
        #include "libmesh/rb_theta.h"
        #include "libmesh/rb_theta_expansion.h"
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
        using libMesh::FEMContext;
        using libMesh::Number;
        using libMesh::Real;
        using libMesh::RBAssemblyExpansion;
        using libMesh::RBParameters;
        using libMesh::RBTheta;
        using libMesh::RBThetaExpansion;
        
</pre>
</div>
<div class = "comment">
boundary IDs
</div>

<div class ="fragment">
<pre>
        #define BOUNDARY_ID_MIN_Z 0
        #define BOUNDARY_ID_MIN_Y 1
        #define BOUNDARY_ID_MAX_X 2
        #define BOUNDARY_ID_MAX_Y 3
        #define BOUNDARY_ID_MIN_X 4
        #define BOUNDARY_ID_MAX_Z 5
        
        class ElasticityRBConstruction;
        
</pre>
</div>
<div class = "comment">
Kronecker delta function
</div>

<div class ="fragment">
<pre>
        inline Real kronecker_delta(unsigned int i,
                                    unsigned int j);
        
</pre>
</div>
<div class = "comment">
Rank-4 tensor for elasticity
</div>

<div class ="fragment">
<pre>
        Real elasticity_tensor(unsigned int i,
                               unsigned int j,
                               unsigned int k,
                               unsigned int l);
        
        struct ElasticityAssembly : ElemAssembly
        {
        
        	ElasticityAssembly(ElasticityRBConstruction& rb_sys_in)
        	:
          rb_sys(rb_sys_in)
          {}
        
        	/**
        	 * The ElasticityRBConstruction object that will use this assembly.
        	 */
        	ElasticityRBConstruction& rb_sys;
        };
        
        struct ThetaA0 : RBTheta { virtual Number evaluate(const RBParameters& ) { return 1.; } };
        struct AssemblyA0 : ElasticityAssembly
        {
        
          AssemblyA0(ElasticityRBConstruction& rb_sys_in)
          :
          ElasticityAssembly(rb_sys_in)
          {}
        
</pre>
</div>
<div class = "comment">
The interior assembly operator
</div>

<div class ="fragment">
<pre>
                virtual void interior_assembly(FEMContext &c);
        
        };
        
        struct ThetaA1 : RBTheta { virtual Number evaluate(const RBParameters& mu) { return mu.get_value("x_scaling"); } };
        struct AssemblyA1 : ElasticityAssembly
        {
        
          AssemblyA1(ElasticityRBConstruction& rb_sys_in)
          :
          ElasticityAssembly(rb_sys_in)
          {}
        
</pre>
</div>
<div class = "comment">
The interior assembly operator
</div>

<div class ="fragment">
<pre>
                virtual void interior_assembly(FEMContext &c);
        
        };
        
        struct ThetaA2 : RBTheta { virtual Number evaluate(const RBParameters& mu) { return 1./mu.get_value("x_scaling"); } };
        struct AssemblyA2 : ElasticityAssembly
        {
        
          AssemblyA2(ElasticityRBConstruction& rb_sys_in)
          :
          ElasticityAssembly(rb_sys_in)
          {}
        
</pre>
</div>
<div class = "comment">
The interior assembly operator
</div>

<div class ="fragment">
<pre>
                virtual void interior_assembly(FEMContext &c);
        
        };
        
        struct ThetaF0 : RBTheta { virtual Number evaluate(const RBParameters& mu) { return mu.get_value("load_Fx"); } };
        struct AssemblyF0 : ElasticityAssembly
        {
          AssemblyF0(ElasticityRBConstruction& rb_sys_in)
          :
          ElasticityAssembly(rb_sys_in)
          {}
        
</pre>
</div>
<div class = "comment">
Apply a traction
</div>

<div class ="fragment">
<pre>
          virtual void boundary_assembly(FEMContext &c);
        
        };
        
        struct ThetaF1 : RBTheta { virtual Number evaluate(const RBParameters& mu)   { return mu.get_value("load_Fy"); } };
        struct AssemblyF1 : ElasticityAssembly
        {
          AssemblyF1(ElasticityRBConstruction& rb_sys_in)
          :
          ElasticityAssembly(rb_sys_in)
          {}
        
</pre>
</div>
<div class = "comment">
Apply a traction
</div>

<div class ="fragment">
<pre>
          virtual void boundary_assembly(FEMContext &c);
        };
        
        struct ThetaF2 : RBTheta { virtual Number evaluate(const RBParameters& mu)   { return mu.get_value("load_Fz"); } };
        struct AssemblyF2 : ElasticityAssembly
        {
          AssemblyF2(ElasticityRBConstruction& rb_sys_in)
          :
          ElasticityAssembly(rb_sys_in)
          {}
        
</pre>
</div>
<div class = "comment">
Apply a traction
</div>

<div class ="fragment">
<pre>
          virtual void boundary_assembly(FEMContext &c);
        };
        
        struct InnerProductAssembly : ElasticityAssembly
        {
        
          InnerProductAssembly(ElasticityRBConstruction& rb_sys_in)
          :
          ElasticityAssembly(rb_sys_in)
          {}
        
</pre>
</div>
<div class = "comment">
The interior assembly operator
</div>

<div class ="fragment">
<pre>
                virtual void interior_assembly(FEMContext &c);
        
        };
        
</pre>
</div>
<div class = "comment">
Define an RBThetaExpansion class for this PDE
</div>

<div class ="fragment">
<pre>
        struct ElasticityThetaExpansion : RBThetaExpansion
        {
        
          /**
           * Constructor.
           */
          ElasticityThetaExpansion()
          {
</pre>
</div>
<div class = "comment">
set up the RBThetaExpansion object
</div>

<div class ="fragment">
<pre>
            attach_A_theta(&theta_a_0);
            attach_A_theta(&theta_a_1);
            attach_A_theta(&theta_a_2);
            attach_F_theta(&theta_f_0);
            attach_F_theta(&theta_f_1);
            attach_F_theta(&theta_f_2);
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
          ThetaF0 theta_f_0;
          ThetaF1 theta_f_1;
          ThetaF2 theta_f_2;
        };
        
</pre>
</div>
<div class = "comment">
Define an RBAssemblyExpansion class for this PDE
</div>

<div class ="fragment">
<pre>
        struct ElasticityAssemblyExpansion : RBAssemblyExpansion
        {
        
          /**
           * Constructor.
           */
          ElasticityAssemblyExpansion(ElasticityRBConstruction& rb_sys_in)
          :
          A0_assembly(rb_sys_in),
          A1_assembly(rb_sys_in),
          A2_assembly(rb_sys_in),
          F0_assembly(rb_sys_in),
          F1_assembly(rb_sys_in),
          F2_assembly(rb_sys_in)
          {
</pre>
</div>
<div class = "comment">
And set up the RBAssemblyExpansion object
</div>

<div class ="fragment">
<pre>
            attach_A_assembly(&A0_assembly);
            attach_A_assembly(&A1_assembly);
            attach_A_assembly(&A2_assembly);
            attach_F_assembly(&F0_assembly);
            attach_F_assembly(&F1_assembly);
            attach_F_assembly(&F2_assembly);
          }
        
</pre>
</div>
<div class = "comment">
The ElemAssembly objects
</div>

<div class ="fragment">
<pre>
          AssemblyA0 A0_assembly;
          AssemblyA1 A1_assembly;
          AssemblyA2 A2_assembly;
          AssemblyF0 F0_assembly;
          AssemblyF1 F1_assembly;
          AssemblyF2 F2_assembly;
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
        
</pre>
</div>
<div class = "comment">
local includes
</div>

<div class ="fragment">
<pre>
        #include "assembly.h"
        
</pre>
</div>
<div class = "comment">
rbOOmit includes
</div>

<div class ="fragment">
<pre>
        #include "libmesh/rb_construction.h"
        
</pre>
</div>
<div class = "comment">
libMesh includes
</div>

<div class ="fragment">
<pre>
        #include "libmesh/fe_base.h"
        #include "libmesh/dof_map.h"
        
        using namespace libMesh;
        
        
        class ElasticityRBEvaluation : public RBEvaluation
        {
        public:
        
          /**
           * Constructor. Just set the theta expansion.
           */
          ElasticityRBEvaluation(const Parallel::Communicator& comm)
            : RBEvaluation(comm)
          {
            set_rb_theta_expansion(elasticity_theta_expansion);
          }
        
          /**
           * Return a "dummy" lower bound for the coercivity constant.
           * To do this rigorously we should use the SCM classes.
           */
          virtual Real get_stability_lower_bound() { return 1.; }
        
          /**
           * The object that stores the "theta" expansion of the parameter dependent PDE,
           * i.e. the set of parameter-dependent functions in the affine expansion of the PDE.
           */
          ElasticityThetaExpansion elasticity_theta_expansion;
        };
        
        
        class ElasticityRBConstruction : public RBConstruction
        {
        public:
        
          ElasticityRBConstruction (EquationSystems& es,
                                    const std::string& name_in,
                                    const unsigned int number_in)
          : Parent(es, name_in, number_in),
            elasticity_assembly_expansion(*this),
            ip_assembly(*this)
          {}
        
          /**
           * Destructor.
           */
          virtual ~ElasticityRBConstruction () {}
        
          /**
           * The type of system.
           */
          typedef ElasticityRBConstruction sys_type;
        
          /**
           * The type of the parent.
           */
          typedef RBConstruction Parent;
        
          /**
           * Initialize data structures.
           */
          virtual void init_data()
          {
            u_var = this-&gt;add_variable("u", FIRST);
            v_var = this-&gt;add_variable("v", FIRST);
            w_var = this-&gt;add_variable("w", FIRST);
        
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
Set the Dirichet boundary condition
</div>

<div class ="fragment">
<pre>
            dirichlet_bc-&gt;b.insert(BOUNDARY_ID_MIN_X); // Dirichlet boundary at x=0
            dirichlet_bc-&gt;variables.push_back(u_var);
            dirichlet_bc-&gt;variables.push_back(v_var);
            dirichlet_bc-&gt;variables.push_back(w_var);
        
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
Set the rb_assembly_expansion for this Construction object
</div>

<div class ="fragment">
<pre>
            set_rb_assembly_expansion(elasticity_assembly_expansion);
        
</pre>
</div>
<div class = "comment">
We need to define an inner product matrix for this problem
</div>

<div class ="fragment">
<pre>
            set_inner_product_assembly(ip_assembly);
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
           * Variable numbers.
           */
          unsigned int u_var;
          unsigned int v_var;
          unsigned int w_var;
        
          /**
           * The object that stores the "assembly" expansion of the parameter dependent PDE.
           */
          ElasticityAssemblyExpansion elasticity_assembly_expansion;
        
          /**
           * Object to assemble the inner product matrix
           */
          InnerProductAssembly ip_assembly;
        
          /**
           * The object that defines which degrees of freedom are on a Dirichlet boundary.
           */
          AutoPtr&lt;DirichletBoundary&gt; dirichlet_bc;
        
        };
        
        #endif
</pre>
</div>

<a name="comments"></a> 
<br><br><br> <h1> The source file assembly.C with comments: </h1> 
<div class = "comment">
</div>

<div class ="fragment">
<pre>
        #include "assembly.h"
        #include "rb_classes.h"
        
</pre>
</div>
<div class = "comment">
libMesh includes
</div>

<div class ="fragment">
<pre>
        #include "libmesh/sparse_matrix.h"
        #include "libmesh/numeric_vector.h"
        #include "libmesh/dense_matrix.h"
        #include "libmesh/dense_submatrix.h"
        #include "libmesh/dense_vector.h"
        #include "libmesh/dense_subvector.h"
        #include "libmesh/fe.h"
        #include "libmesh/fe_interface.h"
        #include "libmesh/fe_base.h"
        #include "libmesh/elem_assembly.h"
        #include "libmesh/quadrature_gauss.h"
        #include "libmesh/boundary_info.h"
        
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
Kronecker delta function
</div>

<div class ="fragment">
<pre>
        inline Real kronecker_delta(unsigned int i,
                                    unsigned int j)
        {
          return i == j ? 1. : 0.;
        }
        
        Real elasticity_tensor(unsigned int i,
                               unsigned int j,
                               unsigned int k,
                               unsigned int l)
        {
</pre>
</div>
<div class = "comment">
Define the Poisson ratio and Young's modulus
</div>

<div class ="fragment">
<pre>
          const Real nu = 0.3;
          const Real E  = 1.;
        
</pre>
</div>
<div class = "comment">
Define the Lame constants (lambda_1 and lambda_2) based on nu and E
</div>

<div class ="fragment">
<pre>
          const Real lambda_1 = E * nu / ( (1. + nu) * (1. - 2.*nu) );
          const Real lambda_2 = 0.5 * E / (1. + nu);
        
          return lambda_1 * kronecker_delta(i,j) * kronecker_delta(k,l)
               + lambda_2 * (kronecker_delta(i,k) * kronecker_delta(j,l) + kronecker_delta(i,l) * kronecker_delta(j,k));
        }
        
        void AssemblyA0::interior_assembly(FEMContext &c)
        {
          const unsigned int n_components = rb_sys.n_vars();
        
</pre>
</div>
<div class = "comment">
make sure we have three components
</div>

<div class ="fragment">
<pre>
          libmesh_assert_equal_to (n_components, 3);
        
          const unsigned int u_var = rb_sys.u_var;
          const unsigned int v_var = rb_sys.v_var;
          const unsigned int w_var = rb_sys.w_var;
        
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
Now we will build the affine operator
</div>

<div class ="fragment">
<pre>
          unsigned int n_qpoints = c.element_qrule-&gt;n_points();
        
          std::vector&lt;unsigned int&gt; n_var_dofs(n_components);
          n_var_dofs[u_var] = c.dof_indices_var[u_var].size();
          n_var_dofs[v_var] = c.dof_indices_var[v_var].size();
          n_var_dofs[w_var] = c.dof_indices_var[w_var].size();
        
          for (unsigned int C_i = 0; C_i &lt; n_components; C_i++)
          {
            unsigned int C_j = 0;
            for (unsigned int C_k = 0; C_k &lt; n_components; C_k++)
            {
              for (unsigned int C_l = 1; C_l &lt; n_components; C_l++)
              {
        
                Real C_ijkl = elasticity_tensor(C_i,C_j,C_k,C_l);
                for (unsigned int qp=0; qp&lt;n_qpoints; qp++)
                {
                  for (unsigned int i=0; i&lt;n_var_dofs[C_i]; i++)
                  {
                    for (unsigned int j=0; j&lt;n_var_dofs[C_k]; j++)
                    {
                      (c.get_elem_jacobian(C_i,C_k))(i,j) +=
                        JxW[qp]*(C_ijkl * dphi[i][qp](C_j)*dphi[j][qp](C_l));
                    }
                  }
                }
        
              }
            }
          }
        
          for (unsigned int C_i = 0; C_i &lt; n_components; C_i++)
          {
            for (unsigned int C_j = 1; C_j &lt; n_components; C_j++)
            {
              for (unsigned int C_k = 0; C_k &lt; n_components; C_k++)
              {
                unsigned int C_l = 0;
        
                Real C_ijkl = elasticity_tensor(C_i,C_j,C_k,C_l);
                for (unsigned int qp=0; qp&lt;n_qpoints; qp++)
                {
                  for (unsigned int i=0; i&lt;n_var_dofs[C_i]; i++)
                  {
                    for (unsigned int j=0; j&lt;n_var_dofs[C_k]; j++)
                    {
                      (c.get_elem_jacobian(C_i,C_k))(i,j) +=
                        JxW[qp]*(C_ijkl * dphi[i][qp](C_j)*dphi[j][qp](C_l));
                    }
                  }
                }
        
              }
            }
          }
        
        }
        
        void AssemblyA1::interior_assembly(FEMContext &c)
        {
          const unsigned int n_components = rb_sys.n_vars();
        
</pre>
</div>
<div class = "comment">
make sure we have three components
</div>

<div class ="fragment">
<pre>
          libmesh_assert_equal_to (n_components, 3);
        
          const unsigned int u_var = rb_sys.u_var;
          const unsigned int v_var = rb_sys.v_var;
          const unsigned int w_var = rb_sys.w_var;
        
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
Now we will build the affine operator
</div>

<div class ="fragment">
<pre>
          unsigned int n_qpoints = c.element_qrule-&gt;n_points();
        
          std::vector&lt;unsigned int&gt; n_var_dofs(n_components);
          n_var_dofs[u_var] = c.dof_indices_var[u_var].size();
          n_var_dofs[v_var] = c.dof_indices_var[v_var].size();
          n_var_dofs[w_var] = c.dof_indices_var[w_var].size();
        
          for (unsigned int C_i = 0; C_i &lt; n_components; C_i++)
          {
            for (unsigned int C_j = 1; C_j &lt; n_components; C_j++)
            {
              for (unsigned int C_k = 0; C_k &lt; n_components; C_k++)
              {
                for (unsigned int C_l = 1; C_l &lt; n_components; C_l++)
                {
        
                  Real C_ijkl = elasticity_tensor(C_i,C_j,C_k,C_l);
                  for (unsigned int qp=0; qp&lt;n_qpoints; qp++)
                  {
                    for (unsigned int i=0; i&lt;n_var_dofs[C_i]; i++)
                    {
                      for (unsigned int j=0; j&lt;n_var_dofs[C_k]; j++)
                      {
                        (c.get_elem_jacobian(C_i,C_k))(i,j) +=
                          JxW[qp]*(C_ijkl * dphi[i][qp](C_j)*dphi[j][qp](C_l));
                      }
                    }
                  }
        
                }
              }
            }
          }
        }
        
        void AssemblyA2::interior_assembly(FEMContext &c)
        {
          const unsigned int n_components = rb_sys.n_vars();
        
</pre>
</div>
<div class = "comment">
make sure we have three components
</div>

<div class ="fragment">
<pre>
          libmesh_assert_equal_to (n_components, 3);
        
          const unsigned int u_var = rb_sys.u_var;
          const unsigned int v_var = rb_sys.v_var;
          const unsigned int w_var = rb_sys.w_var;
        
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
Now we will build the affine operator
</div>

<div class ="fragment">
<pre>
          unsigned int n_qpoints = c.element_qrule-&gt;n_points();
        
          std::vector&lt;unsigned int&gt; n_var_dofs(n_components);
          n_var_dofs[u_var] = c.dof_indices_var[u_var].size();
          n_var_dofs[v_var] = c.dof_indices_var[v_var].size();
          n_var_dofs[w_var] = c.dof_indices_var[w_var].size();
        
          for (unsigned int C_i = 0; C_i &lt; n_components; C_i++)
          {
            unsigned int C_j = 0;
        
            for (unsigned int C_k = 0; C_k &lt; n_components; C_k++)
            {
              unsigned int C_l = 0;
        
              Real C_ijkl = elasticity_tensor(C_i,C_j,C_k,C_l);
              for (unsigned int qp=0; qp&lt;n_qpoints; qp++)
              {
                for (unsigned int i=0; i&lt;n_var_dofs[C_i]; i++)
                {
                  for (unsigned int j=0; j&lt;n_var_dofs[C_k]; j++)
                  {
                    (c.get_elem_jacobian(C_i,C_k))(i,j) +=
                      JxW[qp]*(C_ijkl * dphi[i][qp](C_j)*dphi[j][qp](C_l));
                  }
                }
              }
        
            }
          }
        }
        
        void AssemblyF0::boundary_assembly(FEMContext &c)
        {
          if(rb_sys.get_mesh().boundary_info-&gt;has_boundary_id(c.elem, c.side, BOUNDARY_ID_MAX_X) )
          {
            const unsigned int u_var = 0;
        
            const std::vector&lt;Real&gt; &JxW_side =
              c.side_fe_var[u_var]-&gt;get_JxW();
        
            const std::vector&lt;std::vector&lt;Real&gt; &gt;& phi_side =
              c.side_fe_var[u_var]-&gt;get_phi();
        
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
            unsigned int n_qpoints = c.side_qrule-&gt;n_points();
            DenseSubVector&lt;Number&gt;& Fu = c.get_elem_residual(u_var);
        
            for (unsigned int qp=0; qp &lt; n_qpoints; qp++)
              for (unsigned int i=0; i &lt; n_u_dofs; i++)
              {
                Fu(i) += JxW_side[qp] * ( 1. * phi_side[i][qp] );
              }
          }
        }
        
        void AssemblyF1::boundary_assembly(FEMContext &c)
        {
          if(rb_sys.get_mesh().boundary_info-&gt;has_boundary_id(c.elem, c.side, BOUNDARY_ID_MAX_X) )
          {
            const unsigned int u_var = 0;
            const unsigned int v_var = 1;
        
            const std::vector&lt;Real&gt; &JxW_side =
              c.side_fe_var[u_var]-&gt;get_JxW();
        
            const std::vector&lt;std::vector&lt;Real&gt; &gt;& phi_side =
              c.side_fe_var[u_var]-&gt;get_phi();
        
</pre>
</div>
<div class = "comment">
The number of local degrees of freedom in each variable
</div>

<div class ="fragment">
<pre>
            const unsigned int n_v_dofs = c.dof_indices_var[v_var].size();
        
</pre>
</div>
<div class = "comment">
Now we will build the affine operator
</div>

<div class ="fragment">
<pre>
            unsigned int n_qpoints = c.side_qrule-&gt;n_points();
            DenseSubVector&lt;Number&gt;& Fv = c.get_elem_residual(v_var);
        
            for (unsigned int qp=0; qp &lt; n_qpoints; qp++)
              for (unsigned int i=0; i &lt; n_v_dofs; i++)
              {
                Fv(i) += JxW_side[qp] * ( 1. * phi_side[i][qp] );
              }
          }
        }
        
        void AssemblyF2::boundary_assembly(FEMContext &c)
        {
          if(rb_sys.get_mesh().boundary_info-&gt;has_boundary_id(c.elem, c.side, BOUNDARY_ID_MAX_X) )
          {
            const unsigned int u_var = 0;
            const unsigned int w_var = 2;
        
            const std::vector&lt;Real&gt; &JxW_side =
              c.side_fe_var[u_var]-&gt;get_JxW();
        
            const std::vector&lt;std::vector&lt;Real&gt; &gt;& phi_side =
              c.side_fe_var[u_var]-&gt;get_phi();
        
</pre>
</div>
<div class = "comment">
The number of local degrees of freedom in each variable
</div>

<div class ="fragment">
<pre>
            const unsigned int n_w_dofs = c.dof_indices_var[w_var].size();
        
</pre>
</div>
<div class = "comment">
Now we will build the affine operator
</div>

<div class ="fragment">
<pre>
            unsigned int n_qpoints = c.side_qrule-&gt;n_points();
            DenseSubVector&lt;Number&gt;& Fw = c.get_elem_residual(w_var);
        
            for (unsigned int qp=0; qp &lt; n_qpoints; qp++)
              for (unsigned int i=0; i &lt; n_w_dofs; i++)
              {
                Fw(i) += JxW_side[qp] * ( 1. * phi_side[i][qp] );
              }
          }
        }
        
        void InnerProductAssembly::interior_assembly(FEMContext &c)
        {
          const unsigned int u_var = rb_sys.u_var;
          const unsigned int v_var = rb_sys.v_var;
          const unsigned int w_var = rb_sys.w_var;
        
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
          const unsigned int n_v_dofs = c.dof_indices_var[v_var].size();
          const unsigned int n_w_dofs = c.dof_indices_var[w_var].size();
        
</pre>
</div>
<div class = "comment">
Now we will build the affine operator
</div>

<div class ="fragment">
<pre>
          unsigned int n_qpoints = (c.get_element_qrule())-&gt;n_points();
        
          DenseSubMatrix&lt;Number&gt;& Kuu = c.get_elem_jacobian(u_var,u_var);
          DenseSubMatrix&lt;Number&gt;& Kvv = c.get_elem_jacobian(v_var,v_var);
          DenseSubMatrix&lt;Number&gt;& Kww = c.get_elem_jacobian(w_var,w_var);
        
          for (unsigned int qp=0; qp&lt;n_qpoints; qp++)
          {
              for (unsigned int i=0; i&lt;n_u_dofs; i++)
                for (unsigned int j=0; j&lt;n_u_dofs; j++)
                {
                  Kuu(i,j) += JxW[qp]*(dphi[i][qp]*dphi[j][qp]);
                }
        
              for (unsigned int i=0; i&lt;n_v_dofs; i++)
                for (unsigned int j=0; j&lt;n_v_dofs; j++)
                {
                  Kvv(i,j) += JxW[qp]*(dphi[i][qp]*dphi[j][qp]);
                }
        
              for (unsigned int i=0; i&lt;n_w_dofs; i++)
                for (unsigned int j=0; j&lt;n_w_dofs; j++)
                {
                  Kww(i,j) += JxW[qp]*(dphi[i][qp]*dphi[j][qp]);
                }
          }
        }
</pre>
</div>

<a name="comments"></a> 
<br><br><br> <h1> The source file reduced_basis_ex5.C with comments: </h1> 
<div class = "comment">
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
        #include "libmesh/quadrature_gauss.h"
        #include "libmesh/libmesh_logging.h"
        
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
        
        
          Mesh mesh (init.comm(), dim);
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
                ElasticityRBEvaluation rb_eval(mesh.comm());
        
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
          std::vector&lt; std::vector&lt;dof_id_type&gt; &gt; dof_indices_var(system.n_vars());
          std::vector&lt;dof_id_type&gt; stress_dof_indices_var;
        
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
                dof_id_type dof_index = stress_dof_indices_var[0];
        
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
            dof_id_type dof_index = stress_dof_indices_var[0];
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
<br><br><br> <h1> The source file assembly.h without comments: </h1> 
<pre> 
  #ifndef __assembly_h__
  #define __assembly_h__
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/rb_parameters.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/rb_theta.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/rb_theta_expansion.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/rb_assembly_expansion.h&quot;</FONT></B>
  
  using libMesh::ElemAssembly;
  using libMesh::FEMContext;
  using libMesh::Number;
  using libMesh::Real;
  using libMesh::RBAssemblyExpansion;
  using libMesh::RBParameters;
  using libMesh::RBTheta;
  using libMesh::RBThetaExpansion;
  
  #define BOUNDARY_ID_MIN_Z 0
  #define BOUNDARY_ID_MIN_Y 1
  #define BOUNDARY_ID_MAX_X 2
  #define BOUNDARY_ID_MAX_Y 3
  #define BOUNDARY_ID_MIN_X 4
  #define BOUNDARY_ID_MAX_Z 5
  
  <B><FONT COLOR="#228B22">class</FONT></B> ElasticityRBConstruction;
  
  <B><FONT COLOR="#228B22">inline</FONT></B> Real kronecker_delta(<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i,
                              <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j);
  
  Real elasticity_tensor(<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i,
                         <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j,
                         <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> k,
                         <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> l);
  
  <B><FONT COLOR="#228B22">struct</FONT></B> ElasticityAssembly : ElemAssembly
  {
  
  	ElasticityAssembly(ElasticityRBConstruction&amp; rb_sys_in)
  	:
    rb_sys(rb_sys_in)
    {}
  
  	<I><FONT COLOR="#B22222">/**
  	 * The ElasticityRBConstruction object that will use this assembly.
  	 */</FONT></I>
  	ElasticityRBConstruction&amp; rb_sys;
  };
  
  <B><FONT COLOR="#228B22">struct</FONT></B> ThetaA0 : RBTheta { <B><FONT COLOR="#228B22">virtual</FONT></B> Number evaluate(<B><FONT COLOR="#228B22">const</FONT></B> RBParameters&amp; ) { <B><FONT COLOR="#A020F0">return</FONT></B> 1.; } };
  <B><FONT COLOR="#228B22">struct</FONT></B> AssemblyA0 : ElasticityAssembly
  {
  
    AssemblyA0(ElasticityRBConstruction&amp; rb_sys_in)
    :
    ElasticityAssembly(rb_sys_in)
    {}
  
  	<B><FONT COLOR="#228B22">virtual</FONT></B> <B><FONT COLOR="#228B22">void</FONT></B> interior_assembly(FEMContext &amp;c);
  
  };
  
  <B><FONT COLOR="#228B22">struct</FONT></B> ThetaA1 : RBTheta { <B><FONT COLOR="#228B22">virtual</FONT></B> Number evaluate(<B><FONT COLOR="#228B22">const</FONT></B> RBParameters&amp; mu) { <B><FONT COLOR="#A020F0">return</FONT></B> mu.get_value(<B><FONT COLOR="#BC8F8F">&quot;x_scaling&quot;</FONT></B>); } };
  <B><FONT COLOR="#228B22">struct</FONT></B> AssemblyA1 : ElasticityAssembly
  {
  
    AssemblyA1(ElasticityRBConstruction&amp; rb_sys_in)
    :
    ElasticityAssembly(rb_sys_in)
    {}
  
  	<B><FONT COLOR="#228B22">virtual</FONT></B> <B><FONT COLOR="#228B22">void</FONT></B> interior_assembly(FEMContext &amp;c);
  
  };
  
  <B><FONT COLOR="#228B22">struct</FONT></B> ThetaA2 : RBTheta { <B><FONT COLOR="#228B22">virtual</FONT></B> Number evaluate(<B><FONT COLOR="#228B22">const</FONT></B> RBParameters&amp; mu) { <B><FONT COLOR="#A020F0">return</FONT></B> 1./mu.get_value(<B><FONT COLOR="#BC8F8F">&quot;x_scaling&quot;</FONT></B>); } };
  <B><FONT COLOR="#228B22">struct</FONT></B> AssemblyA2 : ElasticityAssembly
  {
  
    AssemblyA2(ElasticityRBConstruction&amp; rb_sys_in)
    :
    ElasticityAssembly(rb_sys_in)
    {}
  
  	<B><FONT COLOR="#228B22">virtual</FONT></B> <B><FONT COLOR="#228B22">void</FONT></B> interior_assembly(FEMContext &amp;c);
  
  };
  
  <B><FONT COLOR="#228B22">struct</FONT></B> ThetaF0 : RBTheta { <B><FONT COLOR="#228B22">virtual</FONT></B> Number evaluate(<B><FONT COLOR="#228B22">const</FONT></B> RBParameters&amp; mu) { <B><FONT COLOR="#A020F0">return</FONT></B> mu.get_value(<B><FONT COLOR="#BC8F8F">&quot;load_Fx&quot;</FONT></B>); } };
  <B><FONT COLOR="#228B22">struct</FONT></B> AssemblyF0 : ElasticityAssembly
  {
    AssemblyF0(ElasticityRBConstruction&amp; rb_sys_in)
    :
    ElasticityAssembly(rb_sys_in)
    {}
  
    <B><FONT COLOR="#228B22">virtual</FONT></B> <B><FONT COLOR="#228B22">void</FONT></B> boundary_assembly(FEMContext &amp;c);
  
  };
  
  <B><FONT COLOR="#228B22">struct</FONT></B> ThetaF1 : RBTheta { <B><FONT COLOR="#228B22">virtual</FONT></B> Number evaluate(<B><FONT COLOR="#228B22">const</FONT></B> RBParameters&amp; mu)   { <B><FONT COLOR="#A020F0">return</FONT></B> mu.get_value(<B><FONT COLOR="#BC8F8F">&quot;load_Fy&quot;</FONT></B>); } };
  <B><FONT COLOR="#228B22">struct</FONT></B> AssemblyF1 : ElasticityAssembly
  {
    AssemblyF1(ElasticityRBConstruction&amp; rb_sys_in)
    :
    ElasticityAssembly(rb_sys_in)
    {}
  
    <B><FONT COLOR="#228B22">virtual</FONT></B> <B><FONT COLOR="#228B22">void</FONT></B> boundary_assembly(FEMContext &amp;c);
  };
  
  <B><FONT COLOR="#228B22">struct</FONT></B> ThetaF2 : RBTheta { <B><FONT COLOR="#228B22">virtual</FONT></B> Number evaluate(<B><FONT COLOR="#228B22">const</FONT></B> RBParameters&amp; mu)   { <B><FONT COLOR="#A020F0">return</FONT></B> mu.get_value(<B><FONT COLOR="#BC8F8F">&quot;load_Fz&quot;</FONT></B>); } };
  <B><FONT COLOR="#228B22">struct</FONT></B> AssemblyF2 : ElasticityAssembly
  {
    AssemblyF2(ElasticityRBConstruction&amp; rb_sys_in)
    :
    ElasticityAssembly(rb_sys_in)
    {}
  
    <B><FONT COLOR="#228B22">virtual</FONT></B> <B><FONT COLOR="#228B22">void</FONT></B> boundary_assembly(FEMContext &amp;c);
  };
  
  <B><FONT COLOR="#228B22">struct</FONT></B> InnerProductAssembly : ElasticityAssembly
  {
  
    InnerProductAssembly(ElasticityRBConstruction&amp; rb_sys_in)
    :
    ElasticityAssembly(rb_sys_in)
    {}
  
  	<B><FONT COLOR="#228B22">virtual</FONT></B> <B><FONT COLOR="#228B22">void</FONT></B> interior_assembly(FEMContext &amp;c);
  
  };
  
  <B><FONT COLOR="#228B22">struct</FONT></B> ElasticityThetaExpansion : RBThetaExpansion
  {
  
    <I><FONT COLOR="#B22222">/**
     * Constructor.
     */</FONT></I>
    ElasticityThetaExpansion()
    {
      attach_A_theta(&amp;theta_a_0);
      attach_A_theta(&amp;theta_a_1);
      attach_A_theta(&amp;theta_a_2);
      attach_F_theta(&amp;theta_f_0);
      attach_F_theta(&amp;theta_f_1);
      attach_F_theta(&amp;theta_f_2);
    }
  
    ThetaA0 theta_a_0;
    ThetaA1 theta_a_1;
    ThetaA2 theta_a_2;
    ThetaF0 theta_f_0;
    ThetaF1 theta_f_1;
    ThetaF2 theta_f_2;
  };
  
  <B><FONT COLOR="#228B22">struct</FONT></B> ElasticityAssemblyExpansion : RBAssemblyExpansion
  {
  
    <I><FONT COLOR="#B22222">/**
     * Constructor.
     */</FONT></I>
    ElasticityAssemblyExpansion(ElasticityRBConstruction&amp; rb_sys_in)
    :
    A0_assembly(rb_sys_in),
    A1_assembly(rb_sys_in),
    A2_assembly(rb_sys_in),
    F0_assembly(rb_sys_in),
    F1_assembly(rb_sys_in),
    F2_assembly(rb_sys_in)
    {
      attach_A_assembly(&amp;A0_assembly);
      attach_A_assembly(&amp;A1_assembly);
      attach_A_assembly(&amp;A2_assembly);
      attach_F_assembly(&amp;F0_assembly);
      attach_F_assembly(&amp;F1_assembly);
      attach_F_assembly(&amp;F2_assembly);
    }
  
    AssemblyA0 A0_assembly;
    AssemblyA1 A1_assembly;
    AssemblyA2 A2_assembly;
    AssemblyF0 F0_assembly;
    AssemblyF1 F1_assembly;
    AssemblyF2 F2_assembly;
  };
  
  #endif
</pre> 
<a name="nocomments"></a> 
<br><br><br> <h1> The source file rb_classes.h without comments: </h1> 
<pre> 
  #ifndef __rb_classes_h__
  #define __rb_classes_h__
  
  #include <B><FONT COLOR="#BC8F8F">&quot;assembly.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/rb_construction.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/fe_base.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/dof_map.h&quot;</FONT></B>
  
  using namespace libMesh;
  
  
  <B><FONT COLOR="#228B22">class</FONT></B> ElasticityRBEvaluation : <B><FONT COLOR="#228B22">public</FONT></B> RBEvaluation
  {
  <B><FONT COLOR="#228B22">public</FONT></B>:
  
    <I><FONT COLOR="#B22222">/**
     * Constructor. Just set the theta expansion.
     */</FONT></I>
    ElasticityRBEvaluation(<B><FONT COLOR="#228B22">const</FONT></B> Parallel::Communicator&amp; comm)
      : RBEvaluation(comm)
    {
      set_rb_theta_expansion(elasticity_theta_expansion);
    }
  
    <I><FONT COLOR="#B22222">/**
     * Return a &quot;dummy&quot; lower bound for the coercivity constant.
     * To do this rigorously we should use the SCM classes.
     */</FONT></I>
    <B><FONT COLOR="#228B22">virtual</FONT></B> Real get_stability_lower_bound() { <B><FONT COLOR="#A020F0">return</FONT></B> 1.; }
  
    <I><FONT COLOR="#B22222">/**
     * The object that stores the &quot;theta&quot; expansion of the parameter dependent PDE,
     * i.e. the set of parameter-dependent functions in the affine expansion of the PDE.
     */</FONT></I>
    ElasticityThetaExpansion elasticity_theta_expansion;
  };
  
  
  <B><FONT COLOR="#228B22">class</FONT></B> ElasticityRBConstruction : <B><FONT COLOR="#228B22">public</FONT></B> RBConstruction
  {
  <B><FONT COLOR="#228B22">public</FONT></B>:
  
    ElasticityRBConstruction (EquationSystems&amp; es,
                              <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp; name_in,
                              <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> number_in)
    : Parent(es, name_in, number_in),
      elasticity_assembly_expansion(*<B><FONT COLOR="#A020F0">this</FONT></B>),
      ip_assembly(*<B><FONT COLOR="#A020F0">this</FONT></B>)
    {}
  
    <I><FONT COLOR="#B22222">/**
     * Destructor.
     */</FONT></I>
    <B><FONT COLOR="#228B22">virtual</FONT></B> ~ElasticityRBConstruction () {}
  
    <I><FONT COLOR="#B22222">/**
     * The type of system.
     */</FONT></I>
    <B><FONT COLOR="#228B22">typedef</FONT></B> ElasticityRBConstruction sys_type;
  
    <I><FONT COLOR="#B22222">/**
     * The type of the parent.
     */</FONT></I>
    <B><FONT COLOR="#228B22">typedef</FONT></B> RBConstruction Parent;
  
    <I><FONT COLOR="#B22222">/**
     * Initialize data structures.
     */</FONT></I>
    <B><FONT COLOR="#228B22">virtual</FONT></B> <B><FONT COLOR="#228B22">void</FONT></B> init_data()
    {
      u_var = <B><FONT COLOR="#A020F0">this</FONT></B>-&gt;add_variable(<B><FONT COLOR="#BC8F8F">&quot;u&quot;</FONT></B>, FIRST);
      v_var = <B><FONT COLOR="#A020F0">this</FONT></B>-&gt;add_variable(<B><FONT COLOR="#BC8F8F">&quot;v&quot;</FONT></B>, FIRST);
      w_var = <B><FONT COLOR="#A020F0">this</FONT></B>-&gt;add_variable(<B><FONT COLOR="#BC8F8F">&quot;w&quot;</FONT></B>, FIRST);
  
      dirichlet_bc = build_zero_dirichlet_boundary_object();
  
      dirichlet_bc-&gt;b.insert(BOUNDARY_ID_MIN_X); <I><FONT COLOR="#B22222">// Dirichlet boundary at x=0
</FONT></I>      dirichlet_bc-&gt;variables.push_back(u_var);
      dirichlet_bc-&gt;variables.push_back(v_var);
      dirichlet_bc-&gt;variables.push_back(w_var);
  
      get_dof_map().add_dirichlet_boundary(*dirichlet_bc);
  
      <B><FONT COLOR="#5F9EA0">Parent</FONT></B>::init_data();
  
      set_rb_assembly_expansion(elasticity_assembly_expansion);
  
      set_inner_product_assembly(ip_assembly);
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
     * Variable numbers.
     */</FONT></I>
    <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> u_var;
    <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> v_var;
    <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> w_var;
  
    <I><FONT COLOR="#B22222">/**
     * The object that stores the &quot;assembly&quot; expansion of the parameter dependent PDE.
     */</FONT></I>
    ElasticityAssemblyExpansion elasticity_assembly_expansion;
  
    <I><FONT COLOR="#B22222">/**
     * Object to assemble the inner product matrix
     */</FONT></I>
    InnerProductAssembly ip_assembly;
  
    <I><FONT COLOR="#B22222">/**
     * The object that defines which degrees of freedom are on a Dirichlet boundary.
     */</FONT></I>
    AutoPtr&lt;DirichletBoundary&gt; dirichlet_bc;
  
  };
  
  #endif
</pre> 
<a name="nocomments"></a> 
<br><br><br> <h1> The source file assembly.C without comments: </h1> 
<pre> 
  #include <B><FONT COLOR="#BC8F8F">&quot;assembly.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;rb_classes.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/sparse_matrix.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/numeric_vector.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/dense_matrix.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/dense_submatrix.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/dense_vector.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/dense_subvector.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/fe.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/fe_interface.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/fe_base.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/elem_assembly.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/quadrature_gauss.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/boundary_info.h&quot;</FONT></B>
  
  using libMesh::ElemAssembly;
  using libMesh::FEInterface;
  using libMesh::FEMContext;
  using libMesh::Number;
  using libMesh::Point;
  using libMesh::RBTheta;
  using libMesh::Real;
  using libMesh::RealGradient;
  
  <B><FONT COLOR="#228B22">inline</FONT></B> Real kronecker_delta(<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i,
                              <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j)
  {
    <B><FONT COLOR="#A020F0">return</FONT></B> i == j ? 1. : 0.;
  }
  
  Real elasticity_tensor(<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i,
                         <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j,
                         <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> k,
                         <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> l)
  {
    <B><FONT COLOR="#228B22">const</FONT></B> Real nu = 0.3;
    <B><FONT COLOR="#228B22">const</FONT></B> Real E  = 1.;
  
    <B><FONT COLOR="#228B22">const</FONT></B> Real lambda_1 = E * nu / ( (1. + nu) * (1. - 2.*nu) );
    <B><FONT COLOR="#228B22">const</FONT></B> Real lambda_2 = 0.5 * E / (1. + nu);
  
    <B><FONT COLOR="#A020F0">return</FONT></B> lambda_1 * kronecker_delta(i,j) * kronecker_delta(k,l)
         + lambda_2 * (kronecker_delta(i,k) * kronecker_delta(j,l) + kronecker_delta(i,l) * kronecker_delta(j,k));
  }
  
  <B><FONT COLOR="#228B22">void</FONT></B> AssemblyA0::interior_assembly(FEMContext &amp;c)
  {
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_components = rb_sys.n_vars();
  
    libmesh_assert_equal_to (n_components, 3);
  
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> u_var = rb_sys.u_var;
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> v_var = rb_sys.v_var;
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> w_var = rb_sys.w_var;
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Real&gt; &amp;JxW =
      c.element_fe_var[u_var]-&gt;get_JxW();
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;RealGradient&gt; &gt;&amp; dphi =
      c.element_fe_var[u_var]-&gt;get_dphi();
  
    <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_qpoints = c.element_qrule-&gt;n_points();
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B>&gt; n_var_dofs(n_components);
    n_var_dofs[u_var] = c.dof_indices_var[u_var].size();
    n_var_dofs[v_var] = c.dof_indices_var[v_var].size();
    n_var_dofs[w_var] = c.dof_indices_var[w_var].size();
  
    <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> C_i = 0; C_i &lt; n_components; C_i++)
    {
      <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> C_j = 0;
      <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> C_k = 0; C_k &lt; n_components; C_k++)
      {
        <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> C_l = 1; C_l &lt; n_components; C_l++)
        {
  
          Real C_ijkl = elasticity_tensor(C_i,C_j,C_k,C_l);
          <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> qp=0; qp&lt;n_qpoints; qp++)
          {
            <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i&lt;n_var_dofs[C_i]; i++)
            {
              <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j=0; j&lt;n_var_dofs[C_k]; j++)
              {
                (c.get_elem_jacobian(C_i,C_k))(i,j) +=
                  JxW[qp]*(C_ijkl * dphi[i][qp](C_j)*dphi[j][qp](C_l));
              }
            }
          }
  
        }
      }
    }
  
    <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> C_i = 0; C_i &lt; n_components; C_i++)
    {
      <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> C_j = 1; C_j &lt; n_components; C_j++)
      {
        <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> C_k = 0; C_k &lt; n_components; C_k++)
        {
          <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> C_l = 0;
  
          Real C_ijkl = elasticity_tensor(C_i,C_j,C_k,C_l);
          <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> qp=0; qp&lt;n_qpoints; qp++)
          {
            <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i&lt;n_var_dofs[C_i]; i++)
            {
              <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j=0; j&lt;n_var_dofs[C_k]; j++)
              {
                (c.get_elem_jacobian(C_i,C_k))(i,j) +=
                  JxW[qp]*(C_ijkl * dphi[i][qp](C_j)*dphi[j][qp](C_l));
              }
            }
          }
  
        }
      }
    }
  
  }
  
  <B><FONT COLOR="#228B22">void</FONT></B> AssemblyA1::interior_assembly(FEMContext &amp;c)
  {
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_components = rb_sys.n_vars();
  
    libmesh_assert_equal_to (n_components, 3);
  
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> u_var = rb_sys.u_var;
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> v_var = rb_sys.v_var;
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> w_var = rb_sys.w_var;
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Real&gt; &amp;JxW =
      c.element_fe_var[u_var]-&gt;get_JxW();
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;RealGradient&gt; &gt;&amp; dphi =
      c.element_fe_var[u_var]-&gt;get_dphi();
  
    <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_qpoints = c.element_qrule-&gt;n_points();
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B>&gt; n_var_dofs(n_components);
    n_var_dofs[u_var] = c.dof_indices_var[u_var].size();
    n_var_dofs[v_var] = c.dof_indices_var[v_var].size();
    n_var_dofs[w_var] = c.dof_indices_var[w_var].size();
  
    <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> C_i = 0; C_i &lt; n_components; C_i++)
    {
      <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> C_j = 1; C_j &lt; n_components; C_j++)
      {
        <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> C_k = 0; C_k &lt; n_components; C_k++)
        {
          <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> C_l = 1; C_l &lt; n_components; C_l++)
          {
  
            Real C_ijkl = elasticity_tensor(C_i,C_j,C_k,C_l);
            <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> qp=0; qp&lt;n_qpoints; qp++)
            {
              <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i&lt;n_var_dofs[C_i]; i++)
              {
                <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j=0; j&lt;n_var_dofs[C_k]; j++)
                {
                  (c.get_elem_jacobian(C_i,C_k))(i,j) +=
                    JxW[qp]*(C_ijkl * dphi[i][qp](C_j)*dphi[j][qp](C_l));
                }
              }
            }
  
          }
        }
      }
    }
  }
  
  <B><FONT COLOR="#228B22">void</FONT></B> AssemblyA2::interior_assembly(FEMContext &amp;c)
  {
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_components = rb_sys.n_vars();
  
    libmesh_assert_equal_to (n_components, 3);
  
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> u_var = rb_sys.u_var;
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> v_var = rb_sys.v_var;
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> w_var = rb_sys.w_var;
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Real&gt; &amp;JxW =
      c.element_fe_var[u_var]-&gt;get_JxW();
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;RealGradient&gt; &gt;&amp; dphi =
      c.element_fe_var[u_var]-&gt;get_dphi();
  
    <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_qpoints = c.element_qrule-&gt;n_points();
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B>&gt; n_var_dofs(n_components);
    n_var_dofs[u_var] = c.dof_indices_var[u_var].size();
    n_var_dofs[v_var] = c.dof_indices_var[v_var].size();
    n_var_dofs[w_var] = c.dof_indices_var[w_var].size();
  
    <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> C_i = 0; C_i &lt; n_components; C_i++)
    {
      <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> C_j = 0;
  
      <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> C_k = 0; C_k &lt; n_components; C_k++)
      {
        <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> C_l = 0;
  
        Real C_ijkl = elasticity_tensor(C_i,C_j,C_k,C_l);
        <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> qp=0; qp&lt;n_qpoints; qp++)
        {
          <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i&lt;n_var_dofs[C_i]; i++)
          {
            <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j=0; j&lt;n_var_dofs[C_k]; j++)
            {
              (c.get_elem_jacobian(C_i,C_k))(i,j) +=
                JxW[qp]*(C_ijkl * dphi[i][qp](C_j)*dphi[j][qp](C_l));
            }
          }
        }
  
      }
    }
  }
  
  <B><FONT COLOR="#228B22">void</FONT></B> AssemblyF0::boundary_assembly(FEMContext &amp;c)
  {
    <B><FONT COLOR="#A020F0">if</FONT></B>(rb_sys.get_mesh().boundary_info-&gt;has_boundary_id(c.elem, c.side, BOUNDARY_ID_MAX_X) )
    {
      <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> u_var = 0;
  
      <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Real&gt; &amp;JxW_side =
        c.side_fe_var[u_var]-&gt;get_JxW();
  
      <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;Real&gt; &gt;&amp; phi_side =
        c.side_fe_var[u_var]-&gt;get_phi();
  
      <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_u_dofs = c.dof_indices_var[u_var].size();
  
      <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_qpoints = c.side_qrule-&gt;n_points();
      DenseSubVector&lt;Number&gt;&amp; Fu = c.get_elem_residual(u_var);
  
      <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> qp=0; qp &lt; n_qpoints; qp++)
        <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i &lt; n_u_dofs; i++)
        {
          Fu(i) += JxW_side[qp] * ( 1. * phi_side[i][qp] );
        }
    }
  }
  
  <B><FONT COLOR="#228B22">void</FONT></B> AssemblyF1::boundary_assembly(FEMContext &amp;c)
  {
    <B><FONT COLOR="#A020F0">if</FONT></B>(rb_sys.get_mesh().boundary_info-&gt;has_boundary_id(c.elem, c.side, BOUNDARY_ID_MAX_X) )
    {
      <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> u_var = 0;
      <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> v_var = 1;
  
      <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Real&gt; &amp;JxW_side =
        c.side_fe_var[u_var]-&gt;get_JxW();
  
      <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;Real&gt; &gt;&amp; phi_side =
        c.side_fe_var[u_var]-&gt;get_phi();
  
      <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_v_dofs = c.dof_indices_var[v_var].size();
  
      <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_qpoints = c.side_qrule-&gt;n_points();
      DenseSubVector&lt;Number&gt;&amp; Fv = c.get_elem_residual(v_var);
  
      <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> qp=0; qp &lt; n_qpoints; qp++)
        <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i &lt; n_v_dofs; i++)
        {
          Fv(i) += JxW_side[qp] * ( 1. * phi_side[i][qp] );
        }
    }
  }
  
  <B><FONT COLOR="#228B22">void</FONT></B> AssemblyF2::boundary_assembly(FEMContext &amp;c)
  {
    <B><FONT COLOR="#A020F0">if</FONT></B>(rb_sys.get_mesh().boundary_info-&gt;has_boundary_id(c.elem, c.side, BOUNDARY_ID_MAX_X) )
    {
      <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> u_var = 0;
      <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> w_var = 2;
  
      <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Real&gt; &amp;JxW_side =
        c.side_fe_var[u_var]-&gt;get_JxW();
  
      <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;Real&gt; &gt;&amp; phi_side =
        c.side_fe_var[u_var]-&gt;get_phi();
  
      <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_w_dofs = c.dof_indices_var[w_var].size();
  
      <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_qpoints = c.side_qrule-&gt;n_points();
      DenseSubVector&lt;Number&gt;&amp; Fw = c.get_elem_residual(w_var);
  
      <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> qp=0; qp &lt; n_qpoints; qp++)
        <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i &lt; n_w_dofs; i++)
        {
          Fw(i) += JxW_side[qp] * ( 1. * phi_side[i][qp] );
        }
    }
  }
  
  <B><FONT COLOR="#228B22">void</FONT></B> InnerProductAssembly::interior_assembly(FEMContext &amp;c)
  {
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> u_var = rb_sys.u_var;
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> v_var = rb_sys.v_var;
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> w_var = rb_sys.w_var;
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Real&gt; &amp;JxW =
      c.element_fe_var[u_var]-&gt;get_JxW();
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;RealGradient&gt; &gt;&amp; dphi =
      c.element_fe_var[u_var]-&gt;get_dphi();
  
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_u_dofs = c.dof_indices_var[u_var].size();
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_v_dofs = c.dof_indices_var[v_var].size();
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_w_dofs = c.dof_indices_var[w_var].size();
  
    <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_qpoints = (c.get_element_qrule())-&gt;n_points();
  
    DenseSubMatrix&lt;Number&gt;&amp; Kuu = c.get_elem_jacobian(u_var,u_var);
    DenseSubMatrix&lt;Number&gt;&amp; Kvv = c.get_elem_jacobian(v_var,v_var);
    DenseSubMatrix&lt;Number&gt;&amp; Kww = c.get_elem_jacobian(w_var,w_var);
  
    <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> qp=0; qp&lt;n_qpoints; qp++)
    {
        <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i&lt;n_u_dofs; i++)
          <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j=0; j&lt;n_u_dofs; j++)
          {
            Kuu(i,j) += JxW[qp]*(dphi[i][qp]*dphi[j][qp]);
          }
  
        <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i&lt;n_v_dofs; i++)
          <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j=0; j&lt;n_v_dofs; j++)
          {
            Kvv(i,j) += JxW[qp]*(dphi[i][qp]*dphi[j][qp]);
          }
  
        <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i&lt;n_w_dofs; i++)
          <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j=0; j&lt;n_w_dofs; j++)
          {
            Kww(i,j) += JxW[qp]*(dphi[i][qp]*dphi[j][qp]);
          }
    }
  }
</pre> 
<a name="nocomments"></a> 
<br><br><br> <h1> The source file reduced_basis_ex5.C without comments: </h1> 
<pre> 
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
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/quadrature_gauss.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/libmesh_logging.h&quot;</FONT></B>
  
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
  
  
    Mesh mesh (init.comm(), dim);
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
  
  	ElasticityRBEvaluation rb_eval(mesh.comm());
  
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
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt; std::vector&lt;dof_id_type&gt; &gt; dof_indices_var(system.n_vars());
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;dof_id_type&gt; stress_dof_indices_var;
  
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
  
          dof_id_type dof_index = stress_dof_indices_var[0];
  
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
      dof_id_type dof_index = stress_dof_indices_var[0];
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
make[4]: Entering directory `/net/spark/workspace/roystgnr/libmesh/git/devel/examples/reduced_basis/reduced_basis_ex5'
***************************************************************
* Running Example reduced_basis_ex5:
*  mpirun -np 4 example-devel -online_mode 0 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc
***************************************************************
 
 Mesh Information:
  mesh_dimension()=3
  spatial_dimension()=3
  n_nodes()=1845
    n_local_nodes()=495
  n_elem()=1280
    n_local_elem()=320
    n_active_elem()=1280
  n_subdomains()=1
  n_partitions()=4
  n_processors()=4
  n_threads()=1
  processor_id()=0

*** Warning, This code is untested, experimental, or likely to see future API changes: ../src/reduced_basis/rb_parametrized.C, line 41, compiled Apr 19 2013 at 11:42:51 ***
 EquationSystems
  n_systems()=2
   System #0, "RBElasticity"
    Type "RBConstruction"
    Variables={ "u" "v" "w" } 
    Finite Element Types="LAGRANGE", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" 
    Approximation Orders="FIRST", "THIRD" 
    n_dofs()=5535
    n_local_dofs()=1485
    n_constrained_dofs()=135
    n_local_constrained_dofs()=135
    n_vectors()=1
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 62.3577
      Average Off-Processor Bandwidth <= 3.17073
      Maximum  On-Processor Bandwidth <= 81
      Maximum Off-Processor Bandwidth <= 27
    DofMap Constraints
      Number of DoF Constraints = 135
      Average DoF Constraint Length= 0
      Number of Node Constraints = 0
   System #1, "StressSystem"
    Type "Explicit"
    Variables={ "sigma_00" "sigma_01" "sigma_02" "sigma_10" "sigma_11" "sigma_12" "sigma_20" "sigma_21" "sigma_22" "vonMises" } 
    Finite Element Types="MONOMIAL", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" 
    Approximation Orders="CONSTANT", "THIRD" 
    n_dofs()=12800
    n_local_dofs()=3200
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
load_Fx: 1.700481e+00
load_Fy: 1.054889e+00
load_Fz: -6.660755e-01
x_scaling: 1.157624e+00

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 1 ----
Performing RB solves on training set
Maximum (relative) error bound is 132.139

Performing truth solve at parameter:
load_Fx: 2.960553e+00
load_Fy: 4.225294e+00
load_Fz: 8.243980e-01
x_scaling: 6.191687e-01

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 2 ----
Performing RB solves on training set
Maximum (relative) error bound is 2.25461

Performing truth solve at parameter:
load_Fx: 3.205047e+00
load_Fy: -1.465196e-01
load_Fz: -7.598128e-03
x_scaling: 1.383947e+00

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 3 ----
Performing RB solves on training set
Maximum (relative) error bound is 0.342533

Performing truth solve at parameter:
load_Fx: -4.790304e+00
load_Fy: -1.557665e+00
load_Fz: -7.968378e-01
x_scaling: 1.817373e+00

Warning: Linear solver may not have converged! Final linear residual = 6.35591e-05, number of iterations = 5000

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 4 ----
Performing RB solves on training set
Maximum (relative) error bound is 0.193761

Performing truth solve at parameter:
load_Fx: -7.656995e-02
load_Fy: -2.298177e+00
load_Fz: 8.872931e-01
x_scaling: 5.158586e-01

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 5 ----
Performing RB solves on training set
Maximum (relative) error bound is 0.0587596

Performing truth solve at parameter:
load_Fx: -4.068503e+00
load_Fy: -3.554122e+00
load_Fz: 5.598076e-01
x_scaling: 1.919026e+00

Warning: Linear solver may not have converged! Final linear residual = 3.7638e-05, number of iterations = 5000

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 6 ----
Performing RB solves on training set
Maximum (relative) error bound is 0.0186201

Performing truth solve at parameter:
load_Fx: -3.473907e+00
load_Fy: 2.320495e+00
load_Fz: 4.228101e-01
x_scaling: 5.041678e-01

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 7 ----
Performing RB solves on training set
Maximum (relative) error bound is 0.0168994

Performing truth solve at parameter:
load_Fx: 1.088789e+00
load_Fy: 1.344968e+00
load_Fz: 6.960229e-02
x_scaling: 1.178860e+00

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 8 ----
Performing RB solves on training set
Maximum (relative) error bound is 0.00124991

Performing truth solve at parameter:
load_Fx: 1.606260e+00
load_Fy: -4.501228e+00
load_Fz: 1.485251e-01
x_scaling: 8.600192e-01

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 9 ----
Performing RB solves on training set
Maximum (relative) error bound is 0.000622091

Specified error tolerance reached.

 -------------------------------------------------------------------------------------------------------------------
| Processor id:   0                                                                                                 |
| Num Processors: 4                                                                                                 |
| Time:           Fri Apr 19 11:56:53 2013                                                                          |
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
| libMesh Performance: Alive time=63.3488, Active time=63.3138                                                      |
 -------------------------------------------------------------------------------------------------------------------
| Event                                 nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                                 w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|-------------------------------------------------------------------------------------------------------------------|
|                                                                                                                   |
|                                                                                                                   |
| DofMap                                                                                                            |
|   add_neighbors_to_send_list()        2         0.0095      0.004727    0.0115      0.005735    0.01     0.02     |
|   build_constraint_matrix()           2240      0.0057      0.000003    0.0057      0.000003    0.01     0.01     |
|   build_sparsity()                    1         0.0058      0.005754    0.0177      0.017723    0.01     0.03     |
|   cnstrn_elem_mat_vec()               2240      0.0115      0.000005    0.0115      0.000005    0.02     0.02     |
|   create_dof_constraints()            2         0.0124      0.006216    0.0435      0.021757    0.02     0.07     |
|   distribute_dofs()                   2         0.0060      0.002978    0.0225      0.011274    0.01     0.04     |
|   dof_indices()                       13248     0.1175      0.000009    0.1175      0.000009    0.19     0.19     |
|   prepare_send_list()                 2         0.0001      0.000031    0.0001      0.000031    0.00     0.00     |
|   reinit()                            2         0.0108      0.005390    0.0108      0.005390    0.02     0.02     |
|                                                                                                                   |
| FE                                                                                                                |
|   compute_shape_functions()           8288      0.1149      0.000014    0.1149      0.000014    0.18     0.18     |
|   init_shape_functions()              3822      0.0549      0.000014    0.0549      0.000014    0.09     0.09     |
|                                                                                                                   |
| FEMap                                                                                                             |
|   compute_affine_map()                8288      0.0372      0.000004    0.0372      0.000004    0.06     0.06     |
|   compute_face_map()                  3808      0.0159      0.000004    0.0159      0.000004    0.03     0.03     |
|   init_face_shape_functions()         14        0.0001      0.000009    0.0001      0.000009    0.00     0.00     |
|   init_reference_to_physical_map()    3822      0.0465      0.000012    0.0465      0.000012    0.07     0.07     |
|                                                                                                                   |
| Mesh                                                                                                              |
|   find_neighbors()                    1         0.0061      0.006109    0.0064      0.006433    0.01     0.01     |
|   renumber_nodes_and_elem()           2         0.0006      0.000322    0.0006      0.000322    0.00     0.00     |
|                                                                                                                   |
| MeshCommunication                                                                                                 |
|   assign_global_indices()             1         0.0316      0.031594    0.0320      0.032022    0.05     0.05     |
|   compute_hilbert_indices()           2         0.0111      0.005565    0.0111      0.005565    0.02     0.02     |
|   find_global_indices()               2         0.0015      0.000750    0.0137      0.006845    0.00     0.02     |
|   parallel_sort()                     2         0.0005      0.000253    0.0007      0.000339    0.00     0.00     |
|                                                                                                                   |
| MeshTools::Generation                                                                                             |
|   build_cube()                        1         0.0018      0.001790    0.0018      0.001790    0.00     0.00     |
|                                                                                                                   |
| MetisPartitioner                                                                                                  |
|   partition()                         1         0.0169      0.016901    0.0237      0.023669    0.03     0.04     |
|                                                                                                                   |
| Parallel                                                                                                          |
|   allgather()                         17        0.0039      0.000230    0.0040      0.000232    0.01     0.01     |
|   barrier()                           1         0.0000      0.000025    0.0000      0.000025    0.00     0.00     |
|   broadcast()                         22        0.0001      0.000005    0.0001      0.000005    0.00     0.00     |
|   max(bool)                           5         0.0000      0.000003    0.0000      0.000003    0.00     0.00     |
|   max(scalar)                         146       0.0009      0.000006    0.0009      0.000006    0.00     0.00     |
|   max(vector)                         34        0.0003      0.000009    0.0009      0.000026    0.00     0.00     |
|   maxloc(scalar)                      10        0.0822      0.008216    0.0822      0.008216    0.13     0.13     |
|   min(bool)                           168       0.0009      0.000006    0.0009      0.000006    0.00     0.00     |
|   min(scalar)                         137       0.0192      0.000140    0.0192      0.000140    0.03     0.03     |
|   min(vector)                         34        0.0003      0.000010    0.0010      0.000028    0.00     0.00     |
|   probe()                             92        0.0012      0.000013    0.0012      0.000013    0.00     0.00     |
|   receive()                           76        0.0007      0.000010    0.0019      0.000025    0.00     0.00     |
|   send()                              64        0.0002      0.000003    0.0002      0.000003    0.00     0.00     |
|   send_receive()                      68        0.0003      0.000005    0.0020      0.000029    0.00     0.00     |
|   sum()                               40        0.0004      0.000010    0.0005      0.000012    0.00     0.00     |
|                                                                                                                   |
| Parallel::Request                                                                                                 |
|   wait()                              64        0.0001      0.000001    0.0001      0.000001    0.00     0.00     |
|                                                                                                                   |
| Partitioner                                                                                                       |
|   set_node_processor_ids()            1         0.0013      0.001330    0.0016      0.001583    0.00     0.00     |
|   set_parent_processor_ids()          1         0.0004      0.000366    0.0004      0.000366    0.00     0.00     |
|                                                                                                                   |
| PetscLinearSolver                                                                                                 |
|   solve()                             39        61.0388     1.565098    61.0388     1.565098    96.41    96.41    |
|                                                                                                                   |
| RBConstruction                                                                                                    |
|   add_scaled_matrix_and_vector()      7         0.7303      0.104322    1.1022      0.157462    1.15     1.74     |
|   clear()                             1         0.0010      0.001011    0.0010      0.001011    0.00     0.00     |
|   compute_Fq_representor_innerprods() 1         0.0083      0.008303    0.7350      0.734995    0.01     1.16     |
|   compute_max_error_bound()           10        0.0202      0.002022    0.2932      0.029319    0.03     0.46     |
|   enrich_RB_space()                   9         0.0156      0.001729    0.0156      0.001729    0.02     0.02     |
|   train_reduced_basis()               1         0.0021      0.002052    61.9863     61.986339   0.00     97.90    |
|   truth_assembly()                    9         0.2726      0.030285    0.2726      0.030285    0.43     0.43     |
|   truth_solve()                       9         0.0049      0.000546    54.3425     6.038050    0.01     85.83    |
|   update_RB_system_matrices()         9         0.0935      0.010394    0.0935      0.010394    0.15     0.15     |
|   update_residual_terms()             9         0.2573      0.028587    6.5045      0.722718    0.41     10.27    |
|                                                                                                                   |
| RBEvaluation                                                                                                      |
|   clear()                             1         0.0001      0.000134    0.0001      0.000134    0.00     0.00     |
|   compute_residual_dual_norm()        2500      0.1636      0.000065    0.1636      0.000065    0.26     0.26     |
|   rb_solve()                          2500      0.0265      0.000011    0.1903      0.000076    0.04     0.30     |
|   resize_data_structures()            1         0.0001      0.000062    0.0001      0.000062    0.00     0.00     |
|   write_offline_data_to_files()       1         0.0180      0.017962    0.0180      0.017962    0.03     0.03     |
|   write_out_basis_functions()         1         0.0005      0.000544    0.0624      0.062383    0.00     0.10     |
|   write_out_vectors()                 1         0.0291      0.029090    0.0618      0.061839    0.05     0.10     |
 -------------------------------------------------------------------------------------------------------------------
| Totals:                               51882     63.3138                                         100.00            |
 -------------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example reduced_basis_ex5:
*  mpirun -np 4 example-devel -online_mode 0 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc
***************************************************************
***************************************************************
* Running Example reduced_basis_ex5:
*  mpirun -np 4 example-devel -online_mode 1 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc
***************************************************************
 
 Mesh Information:
  mesh_dimension()=3
  spatial_dimension()=3
  n_nodes()=1845
    n_local_nodes()=495
  n_elem()=1280
    n_local_elem()=320
    n_active_elem()=1280
  n_subdomains()=1
  n_partitions()=4
  n_processors()=4
  n_threads()=1
  processor_id()=0

*** Warning, This code is untested, experimental, or likely to see future API changes: ../src/reduced_basis/rb_parametrized.C, line 41, compiled Apr 19 2013 at 11:42:51 ***
 EquationSystems
  n_systems()=2
   System #0, "RBElasticity"
    Type "RBConstruction"
    Variables={ "u" "v" "w" } 
    Finite Element Types="LAGRANGE", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" 
    Approximation Orders="FIRST", "THIRD" 
    n_dofs()=5535
    n_local_dofs()=1485
    n_constrained_dofs()=135
    n_local_constrained_dofs()=135
    n_vectors()=1
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 62.3577
      Average Off-Processor Bandwidth <= 3.17073
      Maximum  On-Processor Bandwidth <= 81
      Maximum Off-Processor Bandwidth <= 27
    DofMap Constraints
      Number of DoF Constraints = 135
      Average DoF Constraint Length= 0
      Number of Node Constraints = 0
   System #1, "StressSystem"
    Type "Explicit"
    Variables={ "sigma_00" "sigma_01" "sigma_02" "sigma_10" "sigma_11" "sigma_12" "sigma_20" "sigma_21" "sigma_22" "vonMises" } 
    Finite Element Types="MONOMIAL", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" 
    Approximation Orders="CONSTANT", "THIRD" 
    n_dofs()=12800
    n_local_dofs()=3200
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

load_Fx: 0.000000e+00
load_Fy: 0.000000e+00
load_Fz: -1.000000e+00
x_scaling: 1.300000e+00


 -------------------------------------------------------------------------------------------------------------------
| Processor id:   0                                                                                                 |
| Num Processors: 4                                                                                                 |
| Time:           Fri Apr 19 11:56:54 2013                                                                          |
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
 ----------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=0.804637, Active time=0.780065                                                 |
 ----------------------------------------------------------------------------------------------------------------
| Event                              nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                              w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|----------------------------------------------------------------------------------------------------------------|
|                                                                                                                |
|                                                                                                                |
| DofMap                                                                                                         |
|   add_neighbors_to_send_list()     2         0.0133      0.006638    0.0158      0.007912    1.70     2.03     |
|   build_sparsity()                 1         0.0105      0.010474    0.0185      0.018478    1.34     2.37     |
|   create_dof_constraints()         2         0.0134      0.006722    0.0453      0.022627    1.72     5.80     |
|   distribute_dofs()                2         0.0082      0.004105    0.0236      0.011777    1.05     3.02     |
|   dof_indices()                    12608     0.0673      0.000005    0.0673      0.000005    8.63     8.63     |
|   prepare_send_list()              2         0.0001      0.000044    0.0001      0.000044    0.01     0.01     |
|   reinit()                         2         0.0146      0.007291    0.0146      0.007291    1.87     1.87     |
|                                                                                                                |
| EquationSystems                                                                                                |
|   build_solution_vector()          1         0.0058      0.005813    0.0264      0.026394    0.75     3.38     |
|                                                                                                                |
| ExodusII_IO                                                                                                    |
|   write_nodal_data()               1         0.4649      0.464935    0.4649      0.464935    59.60    59.60    |
|                                                                                                                |
| FE                                                                                                             |
|   compute_shape_functions()        320       0.0016      0.000005    0.0016      0.000005    0.20     0.20     |
|   init_shape_functions()           1         0.0000      0.000038    0.0000      0.000038    0.00     0.00     |
|                                                                                                                |
| FEMap                                                                                                          |
|   compute_affine_map()             320       0.0021      0.000006    0.0021      0.000006    0.27     0.27     |
|   init_reference_to_physical_map() 1         0.0001      0.000056    0.0001      0.000056    0.01     0.01     |
|                                                                                                                |
| Mesh                                                                                                           |
|   find_neighbors()                 1         0.0065      0.006538    0.0067      0.006671    0.84     0.86     |
|   renumber_nodes_and_elem()        2         0.0007      0.000335    0.0007      0.000335    0.09     0.09     |
|                                                                                                                |
| MeshCommunication                                                                                              |
|   assign_global_indices()          1         0.0320      0.032003    0.0325      0.032479    4.10     4.16     |
|   compute_hilbert_indices()        2         0.0111      0.005569    0.0111      0.005569    1.43     1.43     |
|   find_global_indices()            2         0.0016      0.000787    0.0138      0.006896    0.20     1.77     |
|   parallel_sort()                  2         0.0007      0.000325    0.0007      0.000371    0.08     0.10     |
|                                                                                                                |
| MeshOutput                                                                                                     |
|   write_equation_systems()         1         0.0001      0.000077    0.4915      0.491478    0.01     63.00    |
|                                                                                                                |
| MeshTools::Generation                                                                                          |
|   build_cube()                     1         0.0019      0.001874    0.0019      0.001874    0.24     0.24     |
|                                                                                                                |
| MetisPartitioner                                                                                               |
|   partition()                      1         0.0172      0.017232    0.0240      0.023964    2.21     3.07     |
|                                                                                                                |
| Parallel                                                                                                       |
|   allgather()                      17        0.0002      0.000012    0.0002      0.000015    0.03     0.03     |
|   barrier()                        1         0.0000      0.000030    0.0000      0.000030    0.00     0.00     |
|   broadcast()                      35        0.0001      0.000003    0.0001      0.000003    0.01     0.01     |
|   max(bool)                        1         0.0000      0.000004    0.0000      0.000004    0.00     0.00     |
|   max(scalar)                      152       0.0008      0.000005    0.0008      0.000005    0.10     0.10     |
|   max(vector)                      34        0.0002      0.000007    0.0006      0.000019    0.03     0.08     |
|   min(bool)                        174       0.0006      0.000004    0.0006      0.000004    0.08     0.08     |
|   min(scalar)                      143       0.0215      0.000150    0.0215      0.000150    2.75     2.75     |
|   min(vector)                      34        0.0003      0.000010    0.0008      0.000023    0.04     0.10     |
|   probe()                          76        0.0006      0.000008    0.0006      0.000008    0.08     0.08     |
|   receive()                        70        0.0003      0.000005    0.0009      0.000013    0.04     0.12     |
|   send()                           70        0.0002      0.000003    0.0002      0.000003    0.03     0.03     |
|   send_receive()                   68        0.0004      0.000006    0.0015      0.000022    0.05     0.20     |
|   sum()                            36        0.0015      0.000042    0.0169      0.000470    0.19     2.17     |
|                                                                                                                |
| Parallel::Request                                                                                              |
|   wait()                           70        0.0006      0.000009    0.0006      0.000009    0.08     0.08     |
|                                                                                                                |
| Partitioner                                                                                                    |
|   set_node_processor_ids()         1         0.0015      0.001520    0.0017      0.001677    0.19     0.21     |
|   set_parent_processor_ids()       1         0.0004      0.000389    0.0004      0.000389    0.05     0.05     |
|                                                                                                                |
| RBConstruction                                                                                                 |
|   clear()                          1         0.0002      0.000227    0.0002      0.000227    0.03     0.03     |
|   load_rb_solution()               1         0.0002      0.000195    0.0002      0.000195    0.02     0.02     |
|                                                                                                                |
| RBEvaluation                                                                                                   |
|   clear()                          1         0.0001      0.000052    0.0001      0.000052    0.01     0.01     |
|   compute_residual_dual_norm()     1         0.0002      0.000239    0.0002      0.000239    0.03     0.03     |
|   rb_solve()                       1         0.0069      0.006908    0.0071      0.007147    0.89     0.92     |
|   read_in_basis_functions()        1         0.0000      0.000021    0.0607      0.060660    0.00     7.78     |
|   read_in_vectors()                1         0.0181      0.018062    0.0606      0.060639    2.32     7.77     |
|   read_offline_data_from_files()   1         0.0006      0.000568    0.0006      0.000637    0.07     0.08     |
|   resize_data_structures()         1         0.0001      0.000069    0.0001      0.000069    0.01     0.01     |
|                                                                                                                |
| main                                                                                                           |
|   compute_stresses()               1         0.0507      0.050662    0.0712      0.071152    6.49     9.12     |
 ----------------------------------------------------------------------------------------------------------------
| Totals:                            14269     0.7801                                          100.00            |
 ----------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example reduced_basis_ex5:
*  mpirun -np 4 example-devel -online_mode 1 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc
***************************************************************
make[4]: Leaving directory `/net/spark/workspace/roystgnr/libmesh/git/devel/examples/reduced_basis/reduced_basis_ex5'
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
