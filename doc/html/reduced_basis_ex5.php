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
          ElasticityRBEvaluation()
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
                                    const std::string& name,
                                    const unsigned int number)
          : Parent(es, name, number),
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
          if(rb_sys.get_mesh().boundary_info-&gt;boundary_id(c.elem, c.side) == BOUNDARY_ID_MAX_X)
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
<h1>Reduced Basis: Example 5 - Reduced Basis Cantilever</h1>
3D cantilever beam using the Reduced Basis Method
David Knezevic, Kyunghoon "K" Lee

<br><br>We consider four parameters in this problem:
x_scaling: scales the length of the cantilever
load_Fx: the traction in the x-direction on the right boundary of the cantilever
load_Fy: the traction in the y-direction on the right boundary of the cantilever
load_Fz: the traction in the z-direction on the right boundary of the cantilever


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
        
        
          Mesh mesh (dim);
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
                ElasticityRBEvaluation rb_eval;
        
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
          const Real x_scaling = mu.get_value("x_scaling");
          
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
    ElasticityRBEvaluation()
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
                              <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp; name,
                              <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> number)
    : Parent(es, name, number),
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
    <B><FONT COLOR="#A020F0">if</FONT></B>(rb_sys.get_mesh().boundary_info-&gt;boundary_id(c.elem, c.side) == BOUNDARY_ID_MAX_X)
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
  
  
    Mesh mesh (dim);
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
  
  	ElasticityRBEvaluation rb_eval;
  
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
    <B><FONT COLOR="#228B22">const</FONT></B> Real x_scaling = mu.get_value(<B><FONT COLOR="#BC8F8F">&quot;x_scaling&quot;</FONT></B>);
    
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
***************************************************************
* Running Example reduced_basis_ex5:
*  mpirun -np 12 example-devel -online_mode 0 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
***************************************************************
 
 Mesh Information:
  mesh_dimension()=3
  spatial_dimension()=3
  n_nodes()=1845
    n_local_nodes()=200
  n_elem()=1280
    n_local_elem()=108
    n_active_elem()=1280
  n_subdomains()=1
  n_partitions()=12
  n_processors()=12
  n_threads()=1
  processor_id()=0

*** Warning, This code is untested, experimental, or likely to see future API changes: src/reduced_basis/rb_parametrized.C, line 41, compiled Jan 31 2013 at 21:51:32 ***
 EquationSystems
  n_systems()=2
   System #0, "RBElasticity"
    Type "RBConstruction"
    Variables={ "u" "v" "w" } 
    Finite Element Types="LAGRANGE", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" 
    Approximation Orders="FIRST", "THIRD" 
    n_dofs()=5535
    n_local_dofs()=600
    n_constrained_dofs()=135
    n_local_constrained_dofs()=0
    n_vectors()=1
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 59.0081
      Average Off-Processor Bandwidth <= 11.2358
      Maximum  On-Processor Bandwidth <= 102
      Maximum Off-Processor Bandwidth <= 63
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
    n_local_dofs()=1080
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

*** Warning, This code is deprecated, and likely to be removed in future library versions! src/mesh/boundary_info.C, line 751, compiled Jan 31 2013 at 21:51:11 ***

---- Performing Greedy basis enrichment ----

---- Basis dimension: 0 ----
Performing RB solves on training set
Maximum (relative) error bound is inf

Performing truth solve at parameter:
load_Fx: 5.842202e-02
load_Fy: -2.353683e+00
load_Fz: 1.403767e+00
x_scaling: 1.942906e+00

Warning: Linear solver may not have converged! Final linear residual = 0.207404, number of iterations = 5000

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 1 ----
Performing RB solves on training set
Maximum (relative) error bound is 40.9035

Performing truth solve at parameter:
load_Fx: -1.336125e-02
load_Fy: 3.751598e+00
load_Fz: 1.558421e+00
x_scaling: 5.711629e-01

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 2 ----
Performing RB solves on training set
Maximum (relative) error bound is 1.28787

Performing truth solve at parameter:
load_Fx: 4.962115e+00
load_Fy: 1.269519e+00
load_Fz: 4.739874e-01
x_scaling: 1.603968e+00

Warning: Linear solver may not have converged! Final linear residual = 0.00145644, number of iterations = 5000

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 3 ----
Performing RB solves on training set
Maximum (relative) error bound is 0.416996

Performing truth solve at parameter:
load_Fx: -3.888662e+00
load_Fy: 1.893687e+00
load_Fz: -6.225015e-01
x_scaling: 5.006494e-01

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 4 ----
Performing RB solves on training set
Maximum (relative) error bound is 0.200457

Performing truth solve at parameter:
load_Fx: 4.521673e+00
load_Fy: 2.626519e+00
load_Fz: -5.749369e-02
x_scaling: 5.037054e-01

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 5 ----
Performing RB solves on training set
Maximum (relative) error bound is 0.0417651

Performing truth solve at parameter:
load_Fx: -5.551802e-01
load_Fy: 2.847258e-01
load_Fz: 2.230697e-03
x_scaling: 1.070787e+00

Warning: Linear solver may not have converged! Final linear residual = 3.06445e-10, number of iterations = 5000

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 6 ----
Performing RB solves on training set
Maximum (relative) error bound is 0.037195

Performing truth solve at parameter:
load_Fx: 1.556717e+00
load_Fy: 5.136455e-01
load_Fz: -4.772965e-02
x_scaling: 8.294728e-01

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 7 ----
Performing RB solves on training set
Maximum (relative) error bound is 0.0148347

Performing truth solve at parameter:
load_Fx: -3.288967e+00
load_Fy: 6.900780e-01
load_Fz: 1.293121e-02
x_scaling: 1.980786e+00

Warning: Linear solver may not have converged! Final linear residual = 0.000365286, number of iterations = 5000

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 8 ----
Performing RB solves on training set
Maximum (relative) error bound is 0.00653984

Performing truth solve at parameter:
load_Fx: -4.516163e+00
load_Fy: -1.082849e-01
load_Fz: -3.301830e-01
x_scaling: 9.178301e-01

Warning: Linear solver may not have converged! Final linear residual = 1.93853e-07, number of iterations = 5000

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 9 ----
Performing RB solves on training set
Maximum (relative) error bound is 0.00306682

Performing truth solve at parameter:
load_Fx: -3.191842e+00
load_Fy: -2.351400e+00
load_Fz: 1.761500e-01
x_scaling: 1.686326e+00

Warning: Linear solver may not have converged! Final linear residual = 0.00051711, number of iterations = 5000

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 10 ----
Performing RB solves on training set
Maximum (relative) error bound is 0.000831384

Specified error tolerance reached.
************************************************************************************************************************
***             WIDEN YOUR WINDOW TO 120 CHARACTERS.  Use 'enscript -r -fCourier9' to print this document            ***
************************************************************************************************************************

---------------------------------------------- PETSc Performance Summary: ----------------------------------------------

/workspace/libmesh/examples/reduced_basis/reduced_basis_ex5/.libs/lt-example-devel on a intel-12. named hbar.ices.utexas.edu with 12 processors, by benkirk Thu Jan 31 22:19:09 2013
Using Petsc Release Version 3.3.0, Patch 2, Fri Jul 13 15:42:00 CDT 2012 

                         Max       Max/Min        Avg      Total 
Time (sec):           2.713e+01      1.00000   2.713e+01
Objects:              6.750e+02      1.00000   6.750e+02
Flops:                1.696e+10      2.64261   1.188e+10  1.426e+11
Flops/sec:            6.253e+08      2.64261   4.381e+08  5.257e+09
MPI Messages:         1.460e+05      1.50000   1.217e+05  1.460e+06
MPI Message Lengths:  1.220e+08      1.69010   8.519e+02  1.244e+09
MPI Reductions:       9.864e+04      1.00000

Flop counting convention: 1 flop = 1 real number operation of type (multiply/divide/add/subtract)
                            e.g., VecAXPY() for real vectors of length N --> 2N flops
                            and VecAXPY() for complex vectors of length N --> 8N flops

Summary of Stages:   ----- Time ------  ----- Flops -----  --- Messages ---  -- Message Lengths --  -- Reductions --
                        Avg     %Total     Avg     %Total   counts   %Total     Avg         %Total   counts   %Total 
 0:      Main Stage: 2.7129e+01 100.0%  1.4260e+11 100.0%  1.460e+06 100.0%  8.519e+02      100.0%  9.864e+04 100.0% 

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

VecDot              1091 1.0 5.4178e-02 3.0 1.31e+06 2.0 0.0e+00 0.0e+00 1.1e+03  0  0  0  0  1   0  0  0  0  1   223
VecMDot            45946 1.0 9.7708e+00 4.9 8.50e+08 2.0 0.0e+00 0.0e+00 4.6e+04 22  5  0  0 47  22  5  0  0 47   803
VecNorm            47543 1.0 2.7077e+00 2.5 5.71e+07 2.0 0.0e+00 0.0e+00 4.8e+04  6  0  0  0 48   6  0  0  0 48   194
VecScale           47570 1.0 7.0647e-02 1.6 2.85e+07 2.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0  3727
VecCopy             1715 1.0 5.7139e-03 2.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecSet             50647 1.0 9.7841e-02 1.6 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecAXPY             3186 1.0 2.4303e-02 2.2 3.82e+06 2.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0  1451
VecMAXPY           47500 1.0 5.2632e-01 1.8 9.06e+08 2.0 0.0e+00 0.0e+00 0.0e+00  2  6  0  0  0   2  6  0  0  0 15883
VecAssemblyBegin     270 1.0 9.3203e-02 2.8 0.00e+00 0.0 9.0e+01 2.0e+03 8.1e+02  0  0  0  0  1   0  0  0  0  1     0
VecAssemblyEnd       270 1.0 2.4939e-04 1.3 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecScatterBegin    48571 1.0 3.1345e-01 1.2 0.00e+00 0.0 1.5e+06 8.5e+02 0.0e+00  1  0100 99  0   1  0100 99  0     0
VecScatterEnd      48571 1.0 1.7529e+00 6.9 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  3  0  0  0  0   3  0  0  0  0     0
VecNormalize       47500 1.0 2.7923e+00 2.4 8.55e+07 2.0 0.0e+00 0.0e+00 4.8e+04  7  1  0  0 48   7  1  0  0 48   282
MatMult            47500 1.0 6.6245e+00 1.6 3.67e+09 2.0 1.4e+06 8.5e+02 0.0e+00 18 23 98 97  0  18 23 98 97  0  5036
MatMultAdd          1028 1.0 6.7295e-02 1.7 8.01e+07 2.0 3.1e+04 8.5e+02 0.0e+00  0  1  2  2  0   0  1  2  2  0 10813
MatSolve           47543 1.0 1.1102e+01 2.6 1.11e+10 3.1 0.0e+00 0.0e+00 0.0e+00 31 63  0  0  0  31 63  0  0  0  8074
MatLUFactorNum        21 1.0 2.0914e-01 3.1 2.67e+08 4.8 0.0e+00 0.0e+00 0.0e+00  0  1  0  0  0   0  1  0  0  0  8785
MatILUFactorSym       21 1.0 7.5102e-01 3.9 0.00e+00 0.0 0.0e+00 0.0e+00 6.3e+01  2  0  0  0  0   2  0  0  0  0     0
MatAssemblyBegin    1180 1.0 2.1333e-01 2.1 0.00e+00 0.0 1.8e+02 3.8e+04 2.4e+03  1  0  0  1  2   1  0  0  1  2     0
MatAssemblyEnd      1180 1.0 7.7324e-02 1.4 0.00e+00 0.0 2.7e+03 2.1e+02 3.7e+02  0  0  0  0  0   0  0  0  0  0     0
MatGetRow          49200 2.0 1.8543e-02 2.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatGetRowIJ           21 1.0 2.6226e-05 2.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatGetOrdering        21 1.0 1.5459e-03 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 8.4e+01  0  0  0  0  0   0  0  0  0  0     0
MatZeroEntries        35 1.0 2.9263e-03 1.8 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatAXPY               41 1.0 1.7149e-01 1.0 0.00e+00 0.0 2.5e+03 2.1e+02 6.6e+02  1  0  0  0  1   1  0  0  0  1     0
KSPGMRESOrthog     45946 1.0 1.0156e+01 4.0 1.70e+09 2.0 0.0e+00 0.0e+00 4.6e+04 23 11  0  0 47  23 11  0  0 47  1545
KSPSetUp              64 1.0 4.5729e-04 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPSolve              43 1.0 2.3098e+01 1.0 1.69e+10 2.6 1.4e+06 8.5e+02 9.4e+04 85 99 98 97 95  85 99 98 97 95  6142
PCSetUp               42 1.0 9.6674e-01 3.6 2.67e+08 4.8 0.0e+00 0.0e+00 1.5e+02  2  1  0  0  0   2  1  0  0  0  1900
PCSetUpOnBlocks       43 1.0 9.6520e-01 3.6 2.67e+08 4.8 0.0e+00 0.0e+00 1.5e+02  2  1  0  0  0   2  1  0  0  0  1903
PCApply            47543 1.0 1.1885e+01 2.3 1.11e+10 3.1 0.0e+00 0.0e+00 0.0e+00 33 63  0  0  0  33 63  0  0  0  7543
------------------------------------------------------------------------------------------------------------------------

Memory usage is given in bytes:

Object Type          Creations   Destructions     Memory  Descendants' Mem.
Reports information only for process 0.

--- Event Stage 0: Main Stage

              Vector   260            260      1278288     0
      Vector Scatter    48             48        49728     0
           Index Set   201            201       321920     0
   IS L to G Mapping     2              2         1128     0
              Matrix   159            159     52491420     0
       Krylov Solver     2              2        19360     0
      Preconditioner     2              2         1784     0
              Viewer     1              0            0     0
========================================================================================================================
Average time to get PetscTime(): 9.53674e-08
Average time for MPI_Barrier(): 4.19617e-06
Average time for zero size MPI_Send(): 1.36693e-05
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
| Time:           Thu Jan 31 22:19:09 2013                                                                             |
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
| libMesh Performance: Alive time=27.3972, Active time=27.0655                                                      |
 -------------------------------------------------------------------------------------------------------------------
| Event                                 nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                                 w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|-------------------------------------------------------------------------------------------------------------------|
|                                                                                                                   |
|                                                                                                                   |
| DofMap                                                                                                            |
|   add_neighbors_to_send_list()        2         0.0592      0.029596    0.1253      0.062630    0.22     0.46     |
|   build_constraint_matrix()           756       0.0089      0.000012    0.0089      0.000012    0.03     0.03     |
|   build_sparsity()                    1         0.0522      0.052167    0.0967      0.096712    0.19     0.36     |
|   cnstrn_elem_mat_vec()               756       0.0015      0.000002    0.0015      0.000002    0.01     0.01     |
|   create_dof_constraints()            2         0.0881      0.044044    0.5015      0.250744    0.33     1.85     |
|   distribute_dofs()                   2         0.0525      0.026234    0.1828      0.091386    0.19     0.68     |
|   dof_indices()                       7244      0.9888      0.000136    0.9888      0.000136    3.65     3.65     |
|   prepare_send_list()                 2         0.0011      0.000552    0.0011      0.000552    0.00     0.00     |
|   reinit()                            2         0.1208      0.060376    0.1208      0.060376    0.45     0.45     |
|                                                                                                                   |
| FE                                                                                                                |
|   compute_shape_functions()           2660      0.3364      0.000126    0.3364      0.000126    1.24     1.24     |
|   init_shape_functions()              1162      0.0657      0.000057    0.0657      0.000057    0.24     0.24     |
|                                                                                                                   |
| FEMap                                                                                                             |
|   compute_affine_map()                2660      0.0469      0.000018    0.0469      0.000018    0.17     0.17     |
|   compute_face_map()                  1148      0.0206      0.000018    0.0206      0.000018    0.08     0.08     |
|   init_face_shape_functions()         14        0.0006      0.000044    0.0006      0.000044    0.00     0.00     |
|   init_reference_to_physical_map()    1162      0.0576      0.000050    0.0576      0.000050    0.21     0.21     |
|                                                                                                                   |
| Mesh                                                                                                              |
|   find_neighbors()                    1         0.0462      0.046168    0.0476      0.047610    0.17     0.18     |
|   renumber_nodes_and_elem()           2         0.0020      0.000992    0.0020      0.000992    0.01     0.01     |
|                                                                                                                   |
| MeshCommunication                                                                                                 |
|   assign_global_indices()             1         0.0627      0.062673    0.0702      0.070243    0.23     0.26     |
|   compute_hilbert_indices()           2         0.0212      0.010601    0.0212      0.010601    0.08     0.08     |
|   find_global_indices()               2         0.0080      0.004023    0.0351      0.017530    0.03     0.13     |
|   parallel_sort()                     2         0.0035      0.001766    0.0045      0.002226    0.01     0.02     |
|                                                                                                                   |
| MeshTools::Generation                                                                                             |
|   build_cube()                        1         0.0103      0.010291    0.0103      0.010291    0.04     0.04     |
|                                                                                                                   |
| MetisPartitioner                                                                                                  |
|   partition()                         1         0.0985      0.098541    0.1154      0.115438    0.36     0.43     |
|                                                                                                                   |
| Parallel                                                                                                          |
|   allgather()                         17        0.0107      0.000630    0.0108      0.000635    0.04     0.04     |
|   barrier()                           1         0.0000      0.000021    0.0000      0.000021    0.00     0.00     |
|   broadcast()                         24        0.0002      0.000008    0.0002      0.000008    0.00     0.00     |
|   max(bool)                           5         0.0000      0.000006    0.0000      0.000006    0.00     0.00     |
|   max(scalar)                         146       0.0014      0.000010    0.0014      0.000010    0.01     0.01     |
|   max(vector)                         34        0.0006      0.000016    0.0015      0.000044    0.00     0.01     |
|   maxloc(scalar)                      11        0.0240      0.002180    0.0240      0.002180    0.09     0.09     |
|   min(bool)                           168       0.0016      0.000009    0.0016      0.000009    0.01     0.01     |
|   min(scalar)                         137       0.0456      0.000333    0.0456      0.000333    0.17     0.17     |
|   min(vector)                         34        0.0007      0.000019    0.0017      0.000051    0.00     0.01     |
|   probe()                             316       0.0052      0.000016    0.0052      0.000016    0.02     0.02     |
|   receive()                           268       0.0020      0.000008    0.0070      0.000026    0.01     0.03     |
|   send()                              224       0.0008      0.000003    0.0008      0.000003    0.00     0.00     |
|   send_receive()                      228       0.0022      0.000009    0.0098      0.000043    0.01     0.04     |
|   sum()                               41        0.0011      0.000027    0.0016      0.000038    0.00     0.01     |
|                                                                                                                   |
| Parallel::Request                                                                                                 |
|   wait()                              224       0.0004      0.000002    0.0004      0.000002    0.00     0.00     |
|                                                                                                                   |
| Partitioner                                                                                                       |
|   set_node_processor_ids()            1         0.0032      0.003170    0.0050      0.005042    0.01     0.02     |
|   set_parent_processor_ids()          1         0.0023      0.002348    0.0023      0.002348    0.01     0.01     |
|                                                                                                                   |
| PetscLinearSolver                                                                                                 |
|   solve()                             43        23.1130     0.537512    23.1130     0.537512    85.40    85.40    |
|                                                                                                                   |
| RBConstruction                                                                                                    |
|   add_scaled_matrix_and_vector()      7         0.8176      0.116796    1.8362      0.262319    3.02     6.78     |
|   clear()                             1         0.0014      0.001390    0.0014      0.001390    0.01     0.01     |
|   compute_Fq_representor_innerprods() 1         0.0143      0.014337    0.2317      0.231681    0.05     0.86     |
|   compute_max_error_bound()           11        0.0257      0.002333    0.3706      0.033693    0.09     1.37     |
|   enrich_RB_space()                   10        0.0121      0.001207    0.0121      0.001207    0.04     0.04     |
|   train_reduced_basis()               1         0.0054      0.005437    23.8562     23.856185   0.02     88.14    |
|   truth_assembly()                    10        0.1316      0.013157    0.1316      0.013157    0.49     0.49     |
|   truth_solve()                       10        0.0048      0.000476    20.6413     2.064125    0.02     76.26    |
|   update_RB_system_matrices()         10        0.0515      0.005145    0.0515      0.005145    0.19     0.19     |
|   update_residual_terms()             10        0.1523      0.015231    2.5431      0.254313    0.56     9.40     |
|                                                                                                                   |
| RBEvaluation                                                                                                      |
|   clear()                             1         0.0002      0.000161    0.0002      0.000161    0.00     0.00     |
|   compute_residual_dual_norm()        924       0.2896      0.000313    0.2896      0.000313    1.07     1.07     |
|   rb_solve()                          924       0.0304      0.000033    0.3203      0.000347    0.11     1.18     |
|   resize_data_structures()            1         0.0005      0.000526    0.0005      0.000526    0.00     0.00     |
|   write_offline_data_to_files()       1         0.0012      0.001183    0.0012      0.001183    0.00     0.00     |
|   write_out_basis_functions()         1         0.0001      0.000060    0.2306      0.230632    0.00     0.85     |
|   write_out_vectors()                 1         0.1589      0.158934    0.2306      0.230571    0.59     0.85     |
 -------------------------------------------------------------------------------------------------------------------
| Totals:                               21432     27.0655                                         100.00            |
 -------------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example reduced_basis_ex5:
*  mpirun -np 12 example-devel -online_mode 0 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
***************************************************************
***************************************************************
* Running Example reduced_basis_ex5:
*  mpirun -np 12 example-devel -online_mode 1 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
***************************************************************
 
 Mesh Information:
  mesh_dimension()=3
  spatial_dimension()=3
  n_nodes()=1845
    n_local_nodes()=200
  n_elem()=1280
    n_local_elem()=108
    n_active_elem()=1280
  n_subdomains()=1
  n_partitions()=12
  n_processors()=12
  n_threads()=1
  processor_id()=0

*** Warning, This code is untested, experimental, or likely to see future API changes: src/reduced_basis/rb_parametrized.C, line 41, compiled Jan 31 2013 at 21:51:32 ***
 EquationSystems
  n_systems()=2
   System #0, "RBElasticity"
    Type "RBConstruction"
    Variables={ "u" "v" "w" } 
    Finite Element Types="LAGRANGE", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" 
    Approximation Orders="FIRST", "THIRD" 
    n_dofs()=5535
    n_local_dofs()=600
    n_constrained_dofs()=135
    n_local_constrained_dofs()=0
    n_vectors()=1
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 59.0081
      Average Off-Processor Bandwidth <= 11.2358
      Maximum  On-Processor Bandwidth <= 102
      Maximum Off-Processor Bandwidth <= 63
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
    n_local_dofs()=1080
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

************************************************************************************************************************
***             WIDEN YOUR WINDOW TO 120 CHARACTERS.  Use 'enscript -r -fCourier9' to print this document            ***
************************************************************************************************************************

---------------------------------------------- PETSc Performance Summary: ----------------------------------------------

/workspace/libmesh/examples/reduced_basis/reduced_basis_ex5/.libs/lt-example-devel on a intel-12. named hbar.ices.utexas.edu with 12 processors, by benkirk Thu Jan 31 22:19:11 2013
Using Petsc Release Version 3.3.0, Patch 2, Fri Jul 13 15:42:00 CDT 2012 

                         Max       Max/Min        Avg      Total 
Time (sec):           1.589e+00      1.00007   1.589e+00
Objects:              3.100e+01      1.00000   3.100e+01
Flops:                1.200e+04      2.00000   9.225e+03  1.107e+05
Flops/sec:            7.550e+03      2.00000   5.804e+03  6.965e+04
MPI Messages:         1.800e+01      1.50000   1.500e+01  1.800e+02
MPI Message Lengths:  1.475e+04      1.50705   8.641e+02  1.555e+05
MPI Reductions:       7.700e+01      1.00000

Flop counting convention: 1 flop = 1 real number operation of type (multiply/divide/add/subtract)
                            e.g., VecAXPY() for real vectors of length N --> 2N flops
                            and VecAXPY() for complex vectors of length N --> 8N flops

Summary of Stages:   ----- Time ------  ----- Flops -----  --- Messages ---  -- Message Lengths --  -- Reductions --
                        Avg     %Total     Avg     %Total   counts   %Total     Avg         %Total   counts   %Total 
 0:      Main Stage: 1.5894e+00 100.0%  1.1070e+05 100.0%  1.800e+02 100.0%  8.641e+02      100.0%  7.600e+01  98.7% 

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

VecCopy                2 1.0 1.5974e-05 2.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecSet                18 1.0 2.3603e-05 1.7 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecAXPY               10 1.0 1.4520e-04 3.0 1.20e+04 2.0 0.0e+00 0.0e+00 0.0e+00  0100  0  0  0   0100  0  0  0   762
VecAssemblyBegin      13 1.0 4.0956e-0292.0 0.00e+00 0.0 0.0e+00 0.0e+00 3.9e+01  1  0  0  0 51   1  0  0  0 51     0
VecAssemblyEnd        13 1.0 4.6015e-05 1.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecScatterBegin        2 1.0 9.9897e-05 1.2 0.00e+00 0.0 6.0e+01 1.7e+03 0.0e+00  0  0 33 67  0   0  0 33 67  0     0
VecScatterEnd          2 1.0 3.2902e-05 2.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatZeroEntries         2 1.0 1.7810e-04 2.4 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
------------------------------------------------------------------------------------------------------------------------

Memory usage is given in bytes:

Object Type          Creations   Destructions     Memory  Descendants' Mem.
Reports information only for process 0.

--- Event Stage 0: Main Stage

              Vector    19             19       129136     0
      Vector Scatter     2              2         2072     0
           Index Set     4              4         6820     0
   IS L to G Mapping     2              2         1128     0
              Matrix     3              3       566696     0
              Viewer     1              0            0     0
========================================================================================================================
Average time to get PetscTime(): 0
Average time for MPI_Barrier(): 4.76837e-06
Average time for zero size MPI_Send(): 1.36693e-05
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
| Time:           Thu Jan 31 22:19:11 2013                                                                             |
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
 ----------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=1.62608, Active time=1.53269                                                   |
 ----------------------------------------------------------------------------------------------------------------
| Event                              nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                              w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|----------------------------------------------------------------------------------------------------------------|
|                                                                                                                |
|                                                                                                                |
| DofMap                                                                                                         |
|   add_neighbors_to_send_list()     2         0.0598      0.029902    0.1256      0.062797    3.90     8.19     |
|   build_sparsity()                 1         0.0535      0.053540    0.1032      0.103195    3.49     6.73     |
|   create_dof_constraints()         2         0.0880      0.043980    0.4880      0.244018    5.74     31.84    |
|   distribute_dofs()                2         0.0524      0.026200    0.1843      0.092156    3.42     12.03    |
|   dof_indices()                    7028      0.6086      0.000087    0.6086      0.000087    39.71    39.71    |
|   prepare_send_list()              2         0.0011      0.000547    0.0011      0.000547    0.07     0.07     |
|   reinit()                         2         0.1187      0.059364    0.1187      0.059364    7.75     7.75     |
|                                                                                                                |
| EquationSystems                                                                                                |
|   build_solution_vector()          1         0.0135      0.013481    0.0728      0.072845    0.88     4.75     |
|                                                                                                                |
| ExodusII_IO                                                                                                    |
|   write_nodal_data()               1         0.0193      0.019272    0.0193      0.019272    1.26     1.26     |
|                                                                                                                |
| FE                                                                                                             |
|   compute_shape_functions()        108       0.0024      0.000022    0.0024      0.000022    0.15     0.15     |
|   init_shape_functions()           1         0.0002      0.000188    0.0002      0.000188    0.01     0.01     |
|                                                                                                                |
| FEMap                                                                                                          |
|   compute_affine_map()             108       0.0021      0.000020    0.0021      0.000020    0.14     0.14     |
|   init_reference_to_physical_map() 1         0.0002      0.000173    0.0002      0.000173    0.01     0.01     |
|                                                                                                                |
| Mesh                                                                                                           |
|   find_neighbors()                 1         0.0454      0.045393    0.0466      0.046619    2.96     3.04     |
|   renumber_nodes_and_elem()        2         0.0020      0.001015    0.0020      0.001015    0.13     0.13     |
|                                                                                                                |
| MeshCommunication                                                                                              |
|   assign_global_indices()          1         0.0628      0.062759    0.0655      0.065537    4.09     4.28     |
|   compute_hilbert_indices()        2         0.0215      0.010769    0.0215      0.010769    1.41     1.41     |
|   find_global_indices()            2         0.0081      0.004035    0.0369      0.018448    0.53     2.41     |
|   parallel_sort()                  2         0.0035      0.001748    0.0059      0.002927    0.23     0.38     |
|                                                                                                                |
| MeshOutput                                                                                                     |
|   write_equation_systems()         1         0.0002      0.000156    0.0924      0.092386    0.01     6.03     |
|                                                                                                                |
| MeshTools::Generation                                                                                          |
|   build_cube()                     1         0.0102      0.010156    0.0102      0.010156    0.66     0.66     |
|                                                                                                                |
| MetisPartitioner                                                                                               |
|   partition()                      1         0.0986      0.098552    0.1173      0.117257    6.43     7.65     |
|                                                                                                                |
| Parallel                                                                                                       |
|   allgather()                      17        0.0069      0.000404    0.0070      0.000409    0.45     0.45     |
|   barrier()                        1         0.0073      0.007314    0.0073      0.007314    0.48     0.48     |
|   broadcast()                      35        0.0004      0.000011    0.0003      0.000010    0.03     0.02     |
|   max(bool)                        1         0.0000      0.000011    0.0000      0.000011    0.00     0.00     |
|   max(scalar)                      152       0.0013      0.000008    0.0013      0.000008    0.08     0.08     |
|   max(vector)                      34        0.0005      0.000015    0.0011      0.000034    0.03     0.07     |
|   min(bool)                        174       0.0011      0.000006    0.0011      0.000006    0.07     0.07     |
|   min(scalar)                      143       0.0459      0.000321    0.0459      0.000321    2.99     2.99     |
|   min(vector)                      34        0.0006      0.000017    0.0015      0.000044    0.04     0.10     |
|   probe()                          268       0.0070      0.000026    0.0070      0.000026    0.46     0.46     |
|   receive()                        246       0.0017      0.000007    0.0087      0.000035    0.11     0.57     |
|   send()                           246       0.0011      0.000004    0.0011      0.000004    0.07     0.07     |
|   send_receive()                   228       0.0022      0.000010    0.0119      0.000052    0.14     0.78     |
|   sum()                            36        0.0043      0.000121    0.0055      0.000154    0.28     0.36     |
|                                                                                                                |
| Parallel::Request                                                                                              |
|   wait()                           246       0.0005      0.000002    0.0005      0.000002    0.03     0.03     |
|                                                                                                                |
| Partitioner                                                                                                    |
|   set_node_processor_ids()         1         0.0032      0.003185    0.0053      0.005347    0.21     0.35     |
|   set_parent_processor_ids()       1         0.0024      0.002371    0.0024      0.002371    0.15     0.15     |
|                                                                                                                |
| RBConstruction                                                                                                 |
|   clear()                          1         0.0006      0.000589    0.0006      0.000589    0.04     0.04     |
|   load_rb_solution()               1         0.0004      0.000424    0.0004      0.000424    0.03     0.03     |
|                                                                                                                |
| RBEvaluation                                                                                                   |
|   clear()                          1         0.0001      0.000149    0.0001      0.000149    0.01     0.01     |
|   compute_residual_dual_norm()     1         0.0007      0.000747    0.0007      0.000747    0.05     0.05     |
|   rb_solve()                       1         0.0094      0.009404    0.0102      0.010151    0.61     0.66     |
|   read_in_basis_functions()        1         0.0001      0.000058    0.1535      0.153460    0.00     10.01    |
|   read_in_vectors()                1         0.0792      0.079153    0.1534      0.153402    5.16     10.01    |
|   read_offline_data_from_files()   1         0.0012      0.001164    0.0017      0.001695    0.08     0.11     |
|   resize_data_structures()         1         0.0005      0.000530    0.0005      0.000530    0.03     0.03     |
|                                                                                                                |
| main                                                                                                           |
|   compute_stresses()               1         0.0823      0.082343    0.1429      0.142934    5.37     9.33     |
 ----------------------------------------------------------------------------------------------------------------
| Totals:                            9145      1.5327                                          100.00            |
 ----------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example reduced_basis_ex5:
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
