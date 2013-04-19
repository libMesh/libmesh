<?php $root=""; ?>
<?php require($root."navigation.php"); ?>
<html>
<head>
  <?php load_style($root); ?>
</head>
 
<body>
 
<?php make_navigation("reduced_basis_ex6",$root)?>
 
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
        #include "libmesh/boundary_info.h"
        
</pre>
</div>
<div class = "comment">
rbOOmit includes
</div>

<div class ="fragment">
<pre>
        #include "libmesh/rb_theta.h"
        #include "libmesh/rb_assembly_expansion.h"
        #include "libmesh/rb_parametrized_function.h"
        #include "libmesh/rb_eim_construction.h"
        #include "libmesh/rb_eim_theta.h"
        
</pre>
</div>
<div class = "comment">
Bring in bits from the libMesh namespace.
Just the bits we're using, since this is a header.
</div>

<div class ="fragment">
<pre>
        using libMesh::boundary_id_type;
        using libMesh::DenseSubMatrix;
        using libMesh::ElemAssembly;
        using libMesh::FEInterface;
        using libMesh::FEMContext;
        using libMesh::Number;
        using libMesh::Point;
        using libMesh::RBAssemblyExpansion;
        using libMesh::RBConstruction;
        using libMesh::RBParameters;
        using libMesh::RBParametrizedFunction;
        using libMesh::RBTheta;
        using libMesh::RBThetaExpansion;
        using libMesh::RBEIMAssembly;
        using libMesh::RBEIMConstruction;
        using libMesh::RBEIMEvaluation;
        using libMesh::RBEIMTheta;
        using libMesh::Real;
        using libMesh::RealGradient;
        
        struct ElemAssemblyWithConstruction : ElemAssembly
        {
          RBConstruction* rb_con;
        };
        
</pre>
</div>
<div class = "comment">
The "x component" of the function we're approximating with EIM
</div>

<div class ="fragment">
<pre>
        struct Gx : public RBParametrizedFunction
        {
          virtual Number evaluate(const RBParameters& mu,
                                  const Point& p)
          {
            Real curvature = mu.get_value("curvature");
            return 1. + curvature*p(0);
          }
        };
        
</pre>
</div>
<div class = "comment">
The "y component" of the function we're approximating with EIM
</div>

<div class ="fragment">
<pre>
        struct Gy : public RBParametrizedFunction
        {
          virtual Number evaluate(const RBParameters& mu,
                                  const Point& p)
          {
            Real curvature = mu.get_value("curvature");
            return 1. + curvature*p(0);
          }
        };
        
</pre>
</div>
<div class = "comment">
The "z component" of the function we're approximating with EIM
</div>

<div class ="fragment">
<pre>
        struct Gz : public RBParametrizedFunction
        {
          virtual Number evaluate(const RBParameters& mu,
                                  const Point& p)
          {
            Real curvature = mu.get_value("curvature");
            return 1./(1. + curvature*p(0));
          }
        };
        
        struct ThetaA0 : RBTheta {
        virtual Number evaluate(const RBParameters& mu)
        {
          return mu.get_value("kappa") * mu.get_value("Bi");
        }
        };
        struct AssemblyA0 : ElemAssemblyWithConstruction
        {
          virtual void boundary_assembly(FEMContext &c)
          {
            const std::vector&lt;boundary_id_type&gt; bc_ids =
              rb_con-&gt;get_mesh().boundary_info-&gt;boundary_ids (c.elem,c.side);
            for (std::vector&lt;boundary_id_type&gt;::const_iterator b =
                 bc_ids.begin(); b != bc_ids.end(); ++b)
              if( *b == 1 || *b == 2 || *b == 3 || *b == 4 )
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
                  unsigned int n_sidepoints = c.side_qrule-&gt;n_points();
        
                  for (unsigned int qp=0; qp != n_sidepoints; qp++)
                    for (unsigned int i=0; i != n_u_dofs; i++)
                      for (unsigned int j=0; j != n_u_dofs; j++)
                        c.get_elem_jacobian()(i,j) += JxW_side[qp] * phi_side[j][qp]*phi_side[i][qp];
        
                  break;
                }
          }
        };
        
        struct ThetaA1 : RBTheta {
        virtual Number evaluate(const RBParameters& mu)
        {
          return mu.get_value("kappa") * mu.get_value("Bi") * mu.get_value("curvature");
        }
        };
        struct AssemblyA1 : ElemAssemblyWithConstruction
        {
          virtual void boundary_assembly(FEMContext &c)
          {
            const std::vector&lt;boundary_id_type&gt; bc_ids =
              rb_con-&gt;get_mesh().boundary_info-&gt;boundary_ids (c.elem,c.side);
            for (std::vector&lt;boundary_id_type&gt;::const_iterator b =
                 bc_ids.begin(); b != bc_ids.end(); ++b)
              if( *b == 1 || *b == 3 ) // y == -0.2, y == 0.2
                {
                  const unsigned int u_var = 0;
        
                  const std::vector&lt;Real&gt; &JxW_side =
                    c.side_fe_var[u_var]-&gt;get_JxW();
        
                  const std::vector&lt;std::vector&lt;Real&gt; &gt;& phi_side =
                    c.side_fe_var[u_var]-&gt;get_phi();
        
                  const std::vector&lt;Point&gt;& xyz =
                    c.side_fe_var[u_var]-&gt;get_xyz();
        
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
                  unsigned int n_sidepoints = c.side_qrule-&gt;n_points();
        
                  for (unsigned int qp=0; qp != n_sidepoints; qp++)
                  {
                    Real x_hat = xyz[qp](0);
        
                    for (unsigned int i=0; i != n_u_dofs; i++)
                      for (unsigned int j=0; j != n_u_dofs; j++)
                        c.get_elem_jacobian()(i,j) += JxW_side[qp] * x_hat * phi_side[j][qp]*phi_side[i][qp];
                  }
        
                  break;
                }
          }
        };
        
        struct ThetaA2 : RBTheta {
        virtual Number evaluate(const RBParameters& mu)
        {
          return 0.2*mu.get_value("kappa") * mu.get_value("Bi") * mu.get_value("curvature");
        }
        };
        struct AssemblyA2 : ElemAssemblyWithConstruction
        {
          virtual void boundary_assembly(FEMContext &c)
          {
            const std::vector&lt;boundary_id_type&gt; bc_ids =
              rb_con-&gt;get_mesh().boundary_info-&gt;boundary_ids (c.elem,c.side);
            for (std::vector&lt;boundary_id_type&gt;::const_iterator b =
                 bc_ids.begin(); b != bc_ids.end(); ++b)
              if( *b == 2 || *b == 4) // x == 0.2, x == -0.2
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
                  unsigned int n_sidepoints = c.side_qrule-&gt;n_points();
        
                  if(*b==2)
                  {
                    for (unsigned int qp=0; qp != n_sidepoints; qp++)
                    {
                      for (unsigned int i=0; i != n_u_dofs; i++)
                        for (unsigned int j=0; j != n_u_dofs; j++)
                      c.get_elem_jacobian()(i,j) += JxW_side[qp] * phi_side[j][qp]*phi_side[i][qp];
                    }
                  }
        
                  if(*b==4)
                  {
                    for (unsigned int qp=0; qp != n_sidepoints; qp++)
                    {
                      for (unsigned int i=0; i != n_u_dofs; i++)
                        for (unsigned int j=0; j != n_u_dofs; j++)
                      c.get_elem_jacobian()(i,j) -= JxW_side[qp] * phi_side[j][qp]*phi_side[i][qp];
                    }
                  }
                }
          }
        };
        
        struct ThetaEIM : RBEIMTheta {
        
        ThetaEIM(RBEIMEvaluation& rb_eim_eval_in, unsigned int index_in)
        :
        RBEIMTheta(rb_eim_eval_in, index_in)
        {}
        
        virtual Number evaluate(const RBParameters& mu)
        {
          return mu.get_value("kappa") * RBEIMTheta::evaluate(mu);
        }
        };
        struct AssemblyEIM : RBEIMAssembly
        {
        
          AssemblyEIM(RBEIMConstruction& rb_eim_con_in,
                      unsigned int basis_function_index_in)
          : RBEIMAssembly(rb_eim_con_in,
                          basis_function_index_in)
          {}
        
          virtual void interior_assembly(FEMContext &c)
          {
</pre>
</div>
<div class = "comment">
PDE variable numbers
</div>

<div class ="fragment">
<pre>
            const unsigned int u_var = 0;
        
</pre>
</div>
<div class = "comment">
EIM variable numbers
</div>

<div class ="fragment">
<pre>
            const unsigned int Gx_var = 0;
            const unsigned int Gy_var = 1;
            const unsigned int Gz_var = 2;
        
            const std::vector&lt;Real&gt; &JxW =
              c.element_fe_var[u_var]-&gt;get_JxW();
        
            const std::vector&lt;std::vector&lt;RealGradient&gt; &gt;& dphi =
              c.element_fe_var[u_var]-&gt;get_dphi();
        
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
        
            std::vector&lt;Number&gt; eim_values_Gx;
            evaluate_basis_function(Gx_var,
                                    *c.elem,
                                    qpoints,
                                    eim_values_Gx);
        
            std::vector&lt;Number&gt; eim_values_Gy;
            evaluate_basis_function(Gy_var,
                                    *c.elem,
                                    qpoints,
                                    eim_values_Gy);
        
            std::vector&lt;Number&gt; eim_values_Gz;
            evaluate_basis_function(Gz_var,
                                    *c.elem,
                                    qpoints,
                                    eim_values_Gz);
        
            for (unsigned int qp=0; qp != n_qpoints; qp++)
            {
              for (unsigned int i=0; i != n_u_dofs; i++)
                for (unsigned int j=0; j != n_u_dofs; j++)
                {
                  c.get_elem_jacobian()(i,j) += JxW[qp] * ( eim_values_Gx[qp]*dphi[i][qp](0)*dphi[j][qp](0) +
                                                            eim_values_Gy[qp]*dphi[i][qp](1)*dphi[j][qp](1) +
                                                            eim_values_Gz[qp]*dphi[i][qp](2)*dphi[j][qp](2) );
                }
            }
          }
        
        };
        
        
        struct ThetaF0 : RBTheta { virtual Number evaluate(const RBParameters&   ) { return 1.; } };
        struct AssemblyF0 : ElemAssembly
        {
        
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
        
        struct ThetaF1 : RBTheta { virtual Number evaluate(const RBParameters& mu) { return mu.get_value("curvature"); } };
        struct AssemblyF1 : ElemAssembly
        {
        
          virtual void interior_assembly(FEMContext &c)
          {
            const unsigned int u_var = 0;
        
            const std::vector&lt;Real&gt; &JxW =
              c.element_fe_var[u_var]-&gt;get_JxW();
        
            const std::vector&lt;std::vector&lt;Real&gt; &gt;& phi =
              c.element_fe_var[u_var]-&gt;get_phi();
        
            const std::vector&lt;Point&gt;& xyz =
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
        
            for (unsigned int qp=0; qp != n_qpoints; qp++)
            {
              Real x_hat = xyz[qp](0);
        
              for (unsigned int i=0; i != n_u_dofs; i++)
                c.get_elem_residual()(i) += JxW[qp] * ( 1.*x_hat*phi[i][qp] );
            }
          }
        };
        
        struct Ex6InnerProduct : ElemAssembly
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
        
        struct Ex6EIMInnerProduct : ElemAssembly
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
            const unsigned int Gx_var = 0;
            const unsigned int Gy_var = 1;
            const unsigned int Gz_var = 2;
        
            const std::vector&lt;Real&gt; &JxW =
              c.element_fe_var[Gx_var]-&gt;get_JxW();
        
            const std::vector&lt;std::vector&lt;Real&gt; &gt;& phi =
              c.element_fe_var[Gx_var]-&gt;get_phi();
        
            const unsigned int n_u_dofs = c.dof_indices_var[Gx_var].size();
        
            unsigned int n_qpoints = (c.get_element_qrule())-&gt;n_points();
        
            DenseSubMatrix&lt;Number&gt;& Kxx = c.get_elem_jacobian(Gx_var,Gx_var);
            DenseSubMatrix&lt;Number&gt;& Kyy = c.get_elem_jacobian(Gy_var,Gy_var);
            DenseSubMatrix&lt;Number&gt;& Kzz = c.get_elem_jacobian(Gz_var,Gz_var);
        
            for (unsigned int qp=0; qp != n_qpoints; qp++)
              for (unsigned int i=0; i != n_u_dofs; i++)
                for (unsigned int j=0; j != n_u_dofs; j++)
                {
                  Kxx(i,j) += JxW[qp] * phi[j][qp]*phi[i][qp];
                  Kyy(i,j) += JxW[qp] * phi[j][qp]*phi[i][qp];
                  Kzz(i,j) += JxW[qp] * phi[j][qp]*phi[i][qp];
                }
          }
        };
        
</pre>
</div>
<div class = "comment">
Define an RBThetaExpansion class for this PDE
The A terms depend on EIM, so we deal with them later
</div>

<div class ="fragment">
<pre>
        struct Ex6ThetaExpansion : RBThetaExpansion
        {
        
          /**
           * Constructor.
           */
          Ex6ThetaExpansion()
          {
            attach_A_theta(&theta_a0);
            attach_A_theta(&theta_a1);
            attach_A_theta(&theta_a2);
            attach_F_theta(&theta_f0); // Attach the rhs theta
            attach_F_theta(&theta_f1);
          }
        
</pre>
</div>
<div class = "comment">
The RBTheta member variables
</div>

<div class ="fragment">
<pre>
          ThetaA0 theta_a0;
          ThetaA1 theta_a1;
          ThetaA2 theta_a2;
          ThetaF0 theta_f0;
          ThetaF1 theta_f1;
        };
        
</pre>
</div>
<div class = "comment">
Define an RBAssemblyExpansion class for this PDE
The A terms depend on EIM, so we deal with them later
</div>

<div class ="fragment">
<pre>
        struct Ex6AssemblyExpansion : RBAssemblyExpansion
        {
        
          /**
           * Constructor.
           */
          Ex6AssemblyExpansion(RBConstruction& rb_con)
          {
</pre>
</div>
<div class = "comment">
Point to the RBConstruction object
</div>

<div class ="fragment">
<pre>
            assembly_a0.rb_con = &rb_con;
            assembly_a1.rb_con = &rb_con;
            assembly_a2.rb_con = &rb_con;
        
            attach_A_assembly(&assembly_a0);
            attach_A_assembly(&assembly_a1);
            attach_A_assembly(&assembly_a2);
            attach_F_assembly(&assembly_f0); // Attach the rhs assembly
            attach_F_assembly(&assembly_f1);
          }
        
</pre>
</div>
<div class = "comment">
The ElemAssembly objects
</div>

<div class ="fragment">
<pre>
          AssemblyA0 assembly_a0;
          AssemblyA1 assembly_a1;
          AssemblyA2 assembly_a2;
          AssemblyF0 assembly_f0;
          AssemblyF1 assembly_f1;
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
        #include "libmesh/rb_eim_evaluation.h"
        #include "assembly.h"
        
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
            attach_parametrized_function(&g_x);
            attach_parametrized_function(&g_y);
            attach_parametrized_function(&g_z);
          }
        
          /**
           * Build a ThetaEIM rather than an RBEIMTheta.
           */
          virtual AutoPtr&lt;RBTheta&gt; build_eim_theta(unsigned int index)
          {
            return AutoPtr&lt;RBTheta&gt;(new ThetaEIM(*this, index));
          }
        
          /**
           * Parametrized functions that we approximate with EIM
           */
          Gx g_x;
          Gy g_y;
          Gz g_z;
        
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
            return AutoPtr&lt;ElemAssembly&gt;(new AssemblyEIM(*this, index));
          }
        
          /**
           * Initialize data structures.
           */
          virtual void init_data()
          {
            Gx_var = this-&gt;add_variable ("x_comp_of_G", FIRST);
            Gy_var = this-&gt;add_variable ("y_comp_of_G", FIRST);
            Gz_var = this-&gt;add_variable ("z_comp_of_G", FIRST);
        
            Parent::init_data();
        
            set_inner_product_assembly(eim_ip);
          }
        
          /**
           * Variable numbers.
           */
          unsigned int Gx_var;
          unsigned int Gy_var;
          unsigned int Gz_var;
        
          /**
           * Inner product assembly object
           */
          Ex6EIMInnerProduct eim_ip;
        
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
          SimpleRBEvaluation(const Parallel::Communicator& comm)
            : RBEvaluation(comm)
          {
            set_rb_theta_expansion(ex6_theta_expansion);
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
          Ex6ThetaExpansion ex6_theta_expansion;
        
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
                                const std::string& name_in,
                                const unsigned int number_in)
          : Parent(es, name_in, number_in),
            ex6_assembly_expansion(*this),
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
            dirichlet_bc-&gt;b.insert(5);
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
            set_rb_assembly_expansion(ex6_assembly_expansion);
        
</pre>
</div>
<div class = "comment">
We need to define an inner product matrix for this problem
</div>

<div class ="fragment">
<pre>
            set_inner_product_assembly(ex6_ip);
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
          Ex6AssemblyExpansion ex6_assembly_expansion;
        
          /**
           * The inner product assembly object
           */
          Ex6InnerProduct ex6_ip;
        
          /**
           * The object that defines which degrees of freedom are on a Dirichlet boundary.
           */
          AutoPtr&lt;DirichletBoundary&gt; dirichlet_bc;
        
        };
        
        #endif
</pre>
</div>

<a name="comments"></a> 
<br><br><br> <h1> The source file reduced_basis_ex6.C with comments: </h1> 
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
        
</pre>
</div>
<div class = "comment">
local includes
</div>

<div class ="fragment">
<pre>
        #include "rb_classes.h"
        #include "eim_classes.h"
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
        void transform_mesh_and_plot(EquationSystems& es, Real curvature, const std::string& filename);
        
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
This is a 3D example
</div>

<div class ="fragment">
<pre>
          libmesh_example_assert(3 == LIBMESH_DIM, "3D support");
        
</pre>
</div>
<div class = "comment">
Parse the input file using GetPot
</div>

<div class ="fragment">
<pre>
          std::string eim_parameters = "eim.in";
          std::string rb_parameters  = "rb.in";
          std::string main_parameters = "reduced_basis_ex6.in";
          GetPot infile(main_parameters);
        
          unsigned int n_elem_xy = infile("n_elem_xy", 1);
          unsigned int n_elem_z  = infile("n_elem_z", 1);
        
</pre>
</div>
<div class = "comment">
Do we write the RB basis functions to disk?
</div>

<div class ="fragment">
<pre>
          bool store_basis_functions = infile("store_basis_functions", true);
        
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
Create a mesh, with dimension to be overridden by build_cube, on
the default MPI communicator.
</div>

<div class ="fragment">
<pre>
          Mesh mesh(init.comm());
        
          MeshTools::Generation::build_cube (mesh,
                                             n_elem_xy, n_elem_xy, n_elem_z,
                                             -0.2, 0.2,
                                             -0.2, 0.2,
                                             0., 3.,
                                             HEX8);
        
</pre>
</div>
<div class = "comment">
Create an equation systems object.
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
          equation_systems.print_info();
          mesh.print_info();
        
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
        
          if(!online_mode) // Perform the Offline stage of the RB method
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
attach the EIM theta objects to the RBEvaluation
</div>

<div class ="fragment">
<pre>
            eim_rb_eval.initialize_eim_theta_objects();
            rb_eval.get_rb_theta_expansion().attach_multiple_A_theta(eim_rb_eval.get_eim_theta_objects());
        
</pre>
</div>
<div class = "comment">
attach the EIM assembly objects to the RBConstruction
</div>

<div class ="fragment">
<pre>
            eim_construction.initialize_eim_assembly_objects();
            rb_construction.get_rb_assembly_expansion().attach_multiple_A_assembly(eim_construction.get_eim_assembly_objects());
        
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
          else // Perform the Online stage of the RB method
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
            rb_eval.get_rb_theta_expansion().attach_multiple_A_theta(eim_rb_eval.get_eim_theta_objects());
        
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
            Real online_curvature = infile("online_curvature", 0.);
            Real online_Bi        = infile("online_Bi", 0.);
            Real online_kappa     = infile("online_kappa", 0.);
            RBParameters online_mu;
            online_mu.set_value("curvature", online_curvature);
            online_mu.set_value("Bi", online_Bi);
            online_mu.set_value("kappa", online_kappa);
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
        
              transform_mesh_and_plot(equation_systems,online_curvature,"RB_sol.e");
            }
          }
        
          return 0;
        }
        
        void transform_mesh_and_plot(EquationSystems& es, Real curvature, const std::string& filename)
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
        
            Real x = (*node)(0);
            Real z = (*node)(2);
        
            (*node)(0) = -1./curvature + (1./curvature + x)*cos(curvature*z);
            (*node)(2) = (1./curvature + x)*sin(curvature*z);
          }
        
        #ifdef LIBMESH_HAVE_EXODUS_API
          ExodusII_IO(mesh).write_equation_systems(filename, es);
        #endif
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
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/boundary_info.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/rb_theta.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/rb_assembly_expansion.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/rb_parametrized_function.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/rb_eim_construction.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/rb_eim_theta.h&quot;</FONT></B>
  
  using libMesh::boundary_id_type;
  using libMesh::DenseSubMatrix;
  using libMesh::ElemAssembly;
  using libMesh::FEInterface;
  using libMesh::FEMContext;
  using libMesh::Number;
  using libMesh::Point;
  using libMesh::RBAssemblyExpansion;
  using libMesh::RBConstruction;
  using libMesh::RBParameters;
  using libMesh::RBParametrizedFunction;
  using libMesh::RBTheta;
  using libMesh::RBThetaExpansion;
  using libMesh::RBEIMAssembly;
  using libMesh::RBEIMConstruction;
  using libMesh::RBEIMEvaluation;
  using libMesh::RBEIMTheta;
  using libMesh::Real;
  using libMesh::RealGradient;
  
  <B><FONT COLOR="#228B22">struct</FONT></B> ElemAssemblyWithConstruction : ElemAssembly
  {
    RBConstruction* rb_con;
  };
  
  <B><FONT COLOR="#228B22">struct</FONT></B> Gx : <B><FONT COLOR="#228B22">public</FONT></B> RBParametrizedFunction
  {
    <B><FONT COLOR="#228B22">virtual</FONT></B> Number evaluate(<B><FONT COLOR="#228B22">const</FONT></B> RBParameters&amp; mu,
                            <B><FONT COLOR="#228B22">const</FONT></B> Point&amp; p)
    {
      Real curvature = mu.get_value(<B><FONT COLOR="#BC8F8F">&quot;curvature&quot;</FONT></B>);
      <B><FONT COLOR="#A020F0">return</FONT></B> 1. + curvature*p(0);
    }
  };
  
  <B><FONT COLOR="#228B22">struct</FONT></B> Gy : <B><FONT COLOR="#228B22">public</FONT></B> RBParametrizedFunction
  {
    <B><FONT COLOR="#228B22">virtual</FONT></B> Number evaluate(<B><FONT COLOR="#228B22">const</FONT></B> RBParameters&amp; mu,
                            <B><FONT COLOR="#228B22">const</FONT></B> Point&amp; p)
    {
      Real curvature = mu.get_value(<B><FONT COLOR="#BC8F8F">&quot;curvature&quot;</FONT></B>);
      <B><FONT COLOR="#A020F0">return</FONT></B> 1. + curvature*p(0);
    }
  };
  
  <B><FONT COLOR="#228B22">struct</FONT></B> Gz : <B><FONT COLOR="#228B22">public</FONT></B> RBParametrizedFunction
  {
    <B><FONT COLOR="#228B22">virtual</FONT></B> Number evaluate(<B><FONT COLOR="#228B22">const</FONT></B> RBParameters&amp; mu,
                            <B><FONT COLOR="#228B22">const</FONT></B> Point&amp; p)
    {
      Real curvature = mu.get_value(<B><FONT COLOR="#BC8F8F">&quot;curvature&quot;</FONT></B>);
      <B><FONT COLOR="#A020F0">return</FONT></B> 1./(1. + curvature*p(0));
    }
  };
  
  <B><FONT COLOR="#228B22">struct</FONT></B> ThetaA0 : RBTheta {
  <B><FONT COLOR="#228B22">virtual</FONT></B> Number evaluate(<B><FONT COLOR="#228B22">const</FONT></B> RBParameters&amp; mu)
  {
    <B><FONT COLOR="#A020F0">return</FONT></B> mu.get_value(<B><FONT COLOR="#BC8F8F">&quot;kappa&quot;</FONT></B>) * mu.get_value(<B><FONT COLOR="#BC8F8F">&quot;Bi&quot;</FONT></B>);
  }
  };
  <B><FONT COLOR="#228B22">struct</FONT></B> AssemblyA0 : ElemAssemblyWithConstruction
  {
    <B><FONT COLOR="#228B22">virtual</FONT></B> <B><FONT COLOR="#228B22">void</FONT></B> boundary_assembly(FEMContext &amp;c)
    {
      <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;boundary_id_type&gt; bc_ids =
        rb_con-&gt;get_mesh().boundary_info-&gt;boundary_ids (c.elem,c.side);
      <B><FONT COLOR="#A020F0">for</FONT></B> (std::vector&lt;boundary_id_type&gt;::const_iterator b =
           bc_ids.begin(); b != bc_ids.end(); ++b)
        <B><FONT COLOR="#A020F0">if</FONT></B>( *b == 1 || *b == 2 || *b == 3 || *b == 4 )
          {
            <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> u_var = 0;
  
            <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Real&gt; &amp;JxW_side =
              c.side_fe_var[u_var]-&gt;get_JxW();
  
            <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;Real&gt; &gt;&amp; phi_side =
              c.side_fe_var[u_var]-&gt;get_phi();
  
            <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_u_dofs = c.dof_indices_var[u_var].size();
  
            <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_sidepoints = c.side_qrule-&gt;n_points();
  
            <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> qp=0; qp != n_sidepoints; qp++)
              <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i != n_u_dofs; i++)
                <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j=0; j != n_u_dofs; j++)
                  c.get_elem_jacobian()(i,j) += JxW_side[qp] * phi_side[j][qp]*phi_side[i][qp];
  
            <B><FONT COLOR="#A020F0">break</FONT></B>;
          }
    }
  };
  
  <B><FONT COLOR="#228B22">struct</FONT></B> ThetaA1 : RBTheta {
  <B><FONT COLOR="#228B22">virtual</FONT></B> Number evaluate(<B><FONT COLOR="#228B22">const</FONT></B> RBParameters&amp; mu)
  {
    <B><FONT COLOR="#A020F0">return</FONT></B> mu.get_value(<B><FONT COLOR="#BC8F8F">&quot;kappa&quot;</FONT></B>) * mu.get_value(<B><FONT COLOR="#BC8F8F">&quot;Bi&quot;</FONT></B>) * mu.get_value(<B><FONT COLOR="#BC8F8F">&quot;curvature&quot;</FONT></B>);
  }
  };
  <B><FONT COLOR="#228B22">struct</FONT></B> AssemblyA1 : ElemAssemblyWithConstruction
  {
    <B><FONT COLOR="#228B22">virtual</FONT></B> <B><FONT COLOR="#228B22">void</FONT></B> boundary_assembly(FEMContext &amp;c)
    {
      <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;boundary_id_type&gt; bc_ids =
        rb_con-&gt;get_mesh().boundary_info-&gt;boundary_ids (c.elem,c.side);
      <B><FONT COLOR="#A020F0">for</FONT></B> (std::vector&lt;boundary_id_type&gt;::const_iterator b =
           bc_ids.begin(); b != bc_ids.end(); ++b)
        <B><FONT COLOR="#A020F0">if</FONT></B>( *b == 1 || *b == 3 ) <I><FONT COLOR="#B22222">// y == -0.2, y == 0.2
</FONT></I>          {
            <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> u_var = 0;
  
            <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Real&gt; &amp;JxW_side =
              c.side_fe_var[u_var]-&gt;get_JxW();
  
            <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;Real&gt; &gt;&amp; phi_side =
              c.side_fe_var[u_var]-&gt;get_phi();
  
            <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Point&gt;&amp; xyz =
              c.side_fe_var[u_var]-&gt;get_xyz();
  
            <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_u_dofs = c.dof_indices_var[u_var].size();
  
            <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_sidepoints = c.side_qrule-&gt;n_points();
  
            <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> qp=0; qp != n_sidepoints; qp++)
            {
              Real x_hat = xyz[qp](0);
  
              <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i != n_u_dofs; i++)
                <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j=0; j != n_u_dofs; j++)
                  c.get_elem_jacobian()(i,j) += JxW_side[qp] * x_hat * phi_side[j][qp]*phi_side[i][qp];
            }
  
            <B><FONT COLOR="#A020F0">break</FONT></B>;
          }
    }
  };
  
  <B><FONT COLOR="#228B22">struct</FONT></B> ThetaA2 : RBTheta {
  <B><FONT COLOR="#228B22">virtual</FONT></B> Number evaluate(<B><FONT COLOR="#228B22">const</FONT></B> RBParameters&amp; mu)
  {
    <B><FONT COLOR="#A020F0">return</FONT></B> 0.2*mu.get_value(<B><FONT COLOR="#BC8F8F">&quot;kappa&quot;</FONT></B>) * mu.get_value(<B><FONT COLOR="#BC8F8F">&quot;Bi&quot;</FONT></B>) * mu.get_value(<B><FONT COLOR="#BC8F8F">&quot;curvature&quot;</FONT></B>);
  }
  };
  <B><FONT COLOR="#228B22">struct</FONT></B> AssemblyA2 : ElemAssemblyWithConstruction
  {
    <B><FONT COLOR="#228B22">virtual</FONT></B> <B><FONT COLOR="#228B22">void</FONT></B> boundary_assembly(FEMContext &amp;c)
    {
      <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;boundary_id_type&gt; bc_ids =
        rb_con-&gt;get_mesh().boundary_info-&gt;boundary_ids (c.elem,c.side);
      <B><FONT COLOR="#A020F0">for</FONT></B> (std::vector&lt;boundary_id_type&gt;::const_iterator b =
           bc_ids.begin(); b != bc_ids.end(); ++b)
        <B><FONT COLOR="#A020F0">if</FONT></B>( *b == 2 || *b == 4) <I><FONT COLOR="#B22222">// x == 0.2, x == -0.2
</FONT></I>          {
            <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> u_var = 0;
  
            <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Real&gt; &amp;JxW_side =
              c.side_fe_var[u_var]-&gt;get_JxW();
  
            <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;Real&gt; &gt;&amp; phi_side =
              c.side_fe_var[u_var]-&gt;get_phi();
  
            <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_u_dofs = c.dof_indices_var[u_var].size();
  
            <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_sidepoints = c.side_qrule-&gt;n_points();
  
            <B><FONT COLOR="#A020F0">if</FONT></B>(*b==2)
            {
              <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> qp=0; qp != n_sidepoints; qp++)
              {
                <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i != n_u_dofs; i++)
                  <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j=0; j != n_u_dofs; j++)
                c.get_elem_jacobian()(i,j) += JxW_side[qp] * phi_side[j][qp]*phi_side[i][qp];
              }
            }
  
            <B><FONT COLOR="#A020F0">if</FONT></B>(*b==4)
            {
              <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> qp=0; qp != n_sidepoints; qp++)
              {
                <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i != n_u_dofs; i++)
                  <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j=0; j != n_u_dofs; j++)
                c.get_elem_jacobian()(i,j) -= JxW_side[qp] * phi_side[j][qp]*phi_side[i][qp];
              }
            }
          }
    }
  };
  
  <B><FONT COLOR="#228B22">struct</FONT></B> ThetaEIM : RBEIMTheta {
  
  ThetaEIM(RBEIMEvaluation&amp; rb_eim_eval_in, <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> index_in)
  :
  RBEIMTheta(rb_eim_eval_in, index_in)
  {}
  
  <B><FONT COLOR="#228B22">virtual</FONT></B> Number evaluate(<B><FONT COLOR="#228B22">const</FONT></B> RBParameters&amp; mu)
  {
    <B><FONT COLOR="#A020F0">return</FONT></B> mu.get_value(<B><FONT COLOR="#BC8F8F">&quot;kappa&quot;</FONT></B>) * RBEIMTheta::evaluate(mu);
  }
  };
  <B><FONT COLOR="#228B22">struct</FONT></B> AssemblyEIM : RBEIMAssembly
  {
  
    AssemblyEIM(RBEIMConstruction&amp; rb_eim_con_in,
                <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> basis_function_index_in)
    : RBEIMAssembly(rb_eim_con_in,
                    basis_function_index_in)
    {}
  
    <B><FONT COLOR="#228B22">virtual</FONT></B> <B><FONT COLOR="#228B22">void</FONT></B> interior_assembly(FEMContext &amp;c)
    {
      <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> u_var = 0;
  
      <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> Gx_var = 0;
      <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> Gy_var = 1;
      <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> Gz_var = 2;
  
      <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Real&gt; &amp;JxW =
        c.element_fe_var[u_var]-&gt;get_JxW();
  
      <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;RealGradient&gt; &gt;&amp; dphi =
        c.element_fe_var[u_var]-&gt;get_dphi();
  
      <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Point&gt;&amp; qpoints =
        c.element_fe_var[u_var]-&gt;get_xyz();
  
      <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_u_dofs = c.dof_indices_var[u_var].size();
  
      <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_qpoints = (c.get_element_qrule())-&gt;n_points();
  
      <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;Number&gt; eim_values_Gx;
      evaluate_basis_function(Gx_var,
                              *c.elem,
                              qpoints,
                              eim_values_Gx);
  
      <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;Number&gt; eim_values_Gy;
      evaluate_basis_function(Gy_var,
                              *c.elem,
                              qpoints,
                              eim_values_Gy);
  
      <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;Number&gt; eim_values_Gz;
      evaluate_basis_function(Gz_var,
                              *c.elem,
                              qpoints,
                              eim_values_Gz);
  
      <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> qp=0; qp != n_qpoints; qp++)
      {
        <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i != n_u_dofs; i++)
          <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j=0; j != n_u_dofs; j++)
          {
            c.get_elem_jacobian()(i,j) += JxW[qp] * ( eim_values_Gx[qp]*dphi[i][qp](0)*dphi[j][qp](0) +
                                                      eim_values_Gy[qp]*dphi[i][qp](1)*dphi[j][qp](1) +
                                                      eim_values_Gz[qp]*dphi[i][qp](2)*dphi[j][qp](2) );
          }
      }
    }
  
  };
  
  
  <B><FONT COLOR="#228B22">struct</FONT></B> ThetaF0 : RBTheta { <B><FONT COLOR="#228B22">virtual</FONT></B> Number evaluate(<B><FONT COLOR="#228B22">const</FONT></B> RBParameters&amp;   ) { <B><FONT COLOR="#A020F0">return</FONT></B> 1.; } };
  <B><FONT COLOR="#228B22">struct</FONT></B> AssemblyF0 : ElemAssembly
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
  
  <B><FONT COLOR="#228B22">struct</FONT></B> ThetaF1 : RBTheta { <B><FONT COLOR="#228B22">virtual</FONT></B> Number evaluate(<B><FONT COLOR="#228B22">const</FONT></B> RBParameters&amp; mu) { <B><FONT COLOR="#A020F0">return</FONT></B> mu.get_value(<B><FONT COLOR="#BC8F8F">&quot;curvature&quot;</FONT></B>); } };
  <B><FONT COLOR="#228B22">struct</FONT></B> AssemblyF1 : ElemAssembly
  {
  
    <B><FONT COLOR="#228B22">virtual</FONT></B> <B><FONT COLOR="#228B22">void</FONT></B> interior_assembly(FEMContext &amp;c)
    {
      <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> u_var = 0;
  
      <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Real&gt; &amp;JxW =
        c.element_fe_var[u_var]-&gt;get_JxW();
  
      <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;Real&gt; &gt;&amp; phi =
        c.element_fe_var[u_var]-&gt;get_phi();
  
      <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Point&gt;&amp; xyz =
        c.element_fe_var[u_var]-&gt;get_xyz();
  
      <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_u_dofs = c.dof_indices_var[u_var].size();
  
      <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_qpoints = (c.get_element_qrule())-&gt;n_points();
  
      <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> qp=0; qp != n_qpoints; qp++)
      {
        Real x_hat = xyz[qp](0);
  
        <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i != n_u_dofs; i++)
          c.get_elem_residual()(i) += JxW[qp] * ( 1.*x_hat*phi[i][qp] );
      }
    }
  };
  
  <B><FONT COLOR="#228B22">struct</FONT></B> Ex6InnerProduct : ElemAssembly
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
  
  <B><FONT COLOR="#228B22">struct</FONT></B> Ex6EIMInnerProduct : ElemAssembly
  {
  
    <B><FONT COLOR="#228B22">virtual</FONT></B> <B><FONT COLOR="#228B22">void</FONT></B> interior_assembly(FEMContext &amp;c)
    {
      <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> Gx_var = 0;
      <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> Gy_var = 1;
      <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> Gz_var = 2;
  
      <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Real&gt; &amp;JxW =
        c.element_fe_var[Gx_var]-&gt;get_JxW();
  
      <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;Real&gt; &gt;&amp; phi =
        c.element_fe_var[Gx_var]-&gt;get_phi();
  
      <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_u_dofs = c.dof_indices_var[Gx_var].size();
  
      <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_qpoints = (c.get_element_qrule())-&gt;n_points();
  
      DenseSubMatrix&lt;Number&gt;&amp; Kxx = c.get_elem_jacobian(Gx_var,Gx_var);
      DenseSubMatrix&lt;Number&gt;&amp; Kyy = c.get_elem_jacobian(Gy_var,Gy_var);
      DenseSubMatrix&lt;Number&gt;&amp; Kzz = c.get_elem_jacobian(Gz_var,Gz_var);
  
      <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> qp=0; qp != n_qpoints; qp++)
        <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i != n_u_dofs; i++)
          <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j=0; j != n_u_dofs; j++)
          {
            Kxx(i,j) += JxW[qp] * phi[j][qp]*phi[i][qp];
            Kyy(i,j) += JxW[qp] * phi[j][qp]*phi[i][qp];
            Kzz(i,j) += JxW[qp] * phi[j][qp]*phi[i][qp];
          }
    }
  };
  
  <B><FONT COLOR="#228B22">struct</FONT></B> Ex6ThetaExpansion : RBThetaExpansion
  {
  
    <I><FONT COLOR="#B22222">/**
     * Constructor.
     */</FONT></I>
    Ex6ThetaExpansion()
    {
      attach_A_theta(&amp;theta_a0);
      attach_A_theta(&amp;theta_a1);
      attach_A_theta(&amp;theta_a2);
      attach_F_theta(&amp;theta_f0); <I><FONT COLOR="#B22222">// Attach the rhs theta
</FONT></I>      attach_F_theta(&amp;theta_f1);
    }
  
    ThetaA0 theta_a0;
    ThetaA1 theta_a1;
    ThetaA2 theta_a2;
    ThetaF0 theta_f0;
    ThetaF1 theta_f1;
  };
  
  <B><FONT COLOR="#228B22">struct</FONT></B> Ex6AssemblyExpansion : RBAssemblyExpansion
  {
  
    <I><FONT COLOR="#B22222">/**
     * Constructor.
     */</FONT></I>
    Ex6AssemblyExpansion(RBConstruction&amp; rb_con)
    {
      assembly_a0.rb_con = &amp;rb_con;
      assembly_a1.rb_con = &amp;rb_con;
      assembly_a2.rb_con = &amp;rb_con;
  
      attach_A_assembly(&amp;assembly_a0);
      attach_A_assembly(&amp;assembly_a1);
      attach_A_assembly(&amp;assembly_a2);
      attach_F_assembly(&amp;assembly_f0); <I><FONT COLOR="#B22222">// Attach the rhs assembly
</FONT></I>      attach_F_assembly(&amp;assembly_f1);
    }
  
    AssemblyA0 assembly_a0;
    AssemblyA1 assembly_a1;
    AssemblyA2 assembly_a2;
    AssemblyF0 assembly_f0;
    AssemblyF1 assembly_f1;
  };
  
  #endif
</pre> 
<a name="nocomments"></a> 
<br><br><br> <h1> The source file eim_classes.h without comments: </h1> 
<pre> 
  #ifndef __eim_classes_h__
  #define __eim_classes_h__
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/rb_eim_construction.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/rb_eim_evaluation.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;assembly.h&quot;</FONT></B>
  
  <B><FONT COLOR="#228B22">class</FONT></B> SimpleEIMEvaluation : <B><FONT COLOR="#228B22">public</FONT></B> RBEIMEvaluation
  {
  <B><FONT COLOR="#228B22">public</FONT></B>:
  
    SimpleEIMEvaluation(<B><FONT COLOR="#228B22">const</FONT></B> Parallel::Communicator&amp; comm)
      : RBEIMEvaluation(comm)
    {
      attach_parametrized_function(&amp;g_x);
      attach_parametrized_function(&amp;g_y);
      attach_parametrized_function(&amp;g_z);
    }
  
    <I><FONT COLOR="#B22222">/**
     * Build a ThetaEIM rather than an RBEIMTheta.
     */</FONT></I>
    <B><FONT COLOR="#228B22">virtual</FONT></B> AutoPtr&lt;RBTheta&gt; build_eim_theta(<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> index)
    {
      <B><FONT COLOR="#A020F0">return</FONT></B> AutoPtr&lt;RBTheta&gt;(<B><FONT COLOR="#A020F0">new</FONT></B> ThetaEIM(*<B><FONT COLOR="#A020F0">this</FONT></B>, index));
    }
  
    <I><FONT COLOR="#B22222">/**
     * Parametrized functions that we approximate with EIM
     */</FONT></I>
    Gx g_x;
    Gy g_y;
    Gz g_z;
  
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
      <B><FONT COLOR="#A020F0">return</FONT></B> AutoPtr&lt;ElemAssembly&gt;(<B><FONT COLOR="#A020F0">new</FONT></B> AssemblyEIM(*<B><FONT COLOR="#A020F0">this</FONT></B>, index));
    }
  
    <I><FONT COLOR="#B22222">/**
     * Initialize data structures.
     */</FONT></I>
    <B><FONT COLOR="#228B22">virtual</FONT></B> <B><FONT COLOR="#228B22">void</FONT></B> init_data()
    {
      Gx_var = <B><FONT COLOR="#A020F0">this</FONT></B>-&gt;add_variable (<B><FONT COLOR="#BC8F8F">&quot;x_comp_of_G&quot;</FONT></B>, FIRST);
      Gy_var = <B><FONT COLOR="#A020F0">this</FONT></B>-&gt;add_variable (<B><FONT COLOR="#BC8F8F">&quot;y_comp_of_G&quot;</FONT></B>, FIRST);
      Gz_var = <B><FONT COLOR="#A020F0">this</FONT></B>-&gt;add_variable (<B><FONT COLOR="#BC8F8F">&quot;z_comp_of_G&quot;</FONT></B>, FIRST);
  
      <B><FONT COLOR="#5F9EA0">Parent</FONT></B>::init_data();
  
      set_inner_product_assembly(eim_ip);
    }
  
    <I><FONT COLOR="#B22222">/**
     * Variable numbers.
     */</FONT></I>
    <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> Gx_var;
    <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> Gy_var;
    <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> Gz_var;
  
    <I><FONT COLOR="#B22222">/**
     * Inner product assembly object
     */</FONT></I>
    Ex6EIMInnerProduct eim_ip;
  
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
    SimpleRBEvaluation(<B><FONT COLOR="#228B22">const</FONT></B> Parallel::Communicator&amp; comm)
      : RBEvaluation(comm)
    {
      set_rb_theta_expansion(ex6_theta_expansion);
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
    Ex6ThetaExpansion ex6_theta_expansion;
  
  };
  
  <B><FONT COLOR="#228B22">class</FONT></B> SimpleRBConstruction : <B><FONT COLOR="#228B22">public</FONT></B> RBConstruction
  {
  <B><FONT COLOR="#228B22">public</FONT></B>:
  
    SimpleRBConstruction (EquationSystems&amp; es,
                          <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp; name_in,
                          <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> number_in)
    : Parent(es, name_in, number_in),
      ex6_assembly_expansion(*<B><FONT COLOR="#A020F0">this</FONT></B>),
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
      dirichlet_bc-&gt;b.insert(5);
      dirichlet_bc-&gt;variables.push_back(u_var);
  
      get_dof_map().add_dirichlet_boundary(*dirichlet_bc);
  
      <B><FONT COLOR="#5F9EA0">Parent</FONT></B>::init_data();
  
      set_rb_assembly_expansion(ex6_assembly_expansion);
  
      set_inner_product_assembly(ex6_ip);
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
    Ex6AssemblyExpansion ex6_assembly_expansion;
  
    <I><FONT COLOR="#B22222">/**
     * The inner product assembly object
     */</FONT></I>
    Ex6InnerProduct ex6_ip;
  
    <I><FONT COLOR="#B22222">/**
     * The object that defines which degrees of freedom are on a Dirichlet boundary.
     */</FONT></I>
    AutoPtr&lt;DirichletBoundary&gt; dirichlet_bc;
  
  };
  
  #endif
</pre> 
<a name="nocomments"></a> 
<br><br><br> <h1> The source file reduced_basis_ex6.C without comments: </h1> 
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
  
  #include <B><FONT COLOR="#BC8F8F">&quot;rb_classes.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;eim_classes.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;assembly.h&quot;</FONT></B>
  
  using namespace libMesh;
  
  <B><FONT COLOR="#228B22">void</FONT></B> transform_mesh_and_plot(EquationSystems&amp; es, Real curvature, <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp; filename);
  
  <B><FONT COLOR="#228B22">int</FONT></B> main (<B><FONT COLOR="#228B22">int</FONT></B> argc, <B><FONT COLOR="#228B22">char</FONT></B>** argv)
  {
    LibMeshInit init (argc, argv);
  
  #<B><FONT COLOR="#A020F0">if</FONT></B> !defined(LIBMESH_HAVE_XDR)
    libmesh_example_assert(false, <B><FONT COLOR="#BC8F8F">&quot;--enable-xdr&quot;</FONT></B>);
  #elif defined(LIBMESH_DEFAULT_SINGLE_PRECISION)
    libmesh_example_assert(false, <B><FONT COLOR="#BC8F8F">&quot;--disable-singleprecision&quot;</FONT></B>);
  #endif
  
    libmesh_example_assert(3 == LIBMESH_DIM, <B><FONT COLOR="#BC8F8F">&quot;3D support&quot;</FONT></B>);
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::string eim_parameters = <B><FONT COLOR="#BC8F8F">&quot;eim.in&quot;</FONT></B>;
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::string rb_parameters  = <B><FONT COLOR="#BC8F8F">&quot;rb.in&quot;</FONT></B>;
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::string main_parameters = <B><FONT COLOR="#BC8F8F">&quot;reduced_basis_ex6.in&quot;</FONT></B>;
    GetPot infile(main_parameters);
  
    <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_elem_xy = infile(<B><FONT COLOR="#BC8F8F">&quot;n_elem_xy&quot;</FONT></B>, 1);
    <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_elem_z  = infile(<B><FONT COLOR="#BC8F8F">&quot;n_elem_z&quot;</FONT></B>, 1);
  
    <B><FONT COLOR="#228B22">bool</FONT></B> store_basis_functions = infile(<B><FONT COLOR="#BC8F8F">&quot;store_basis_functions&quot;</FONT></B>, true);
  
    GetPot command_line (argc, argv);
    <B><FONT COLOR="#228B22">int</FONT></B> online_mode = 0;
    <B><FONT COLOR="#A020F0">if</FONT></B> ( command_line.search(1, <B><FONT COLOR="#BC8F8F">&quot;-online_mode&quot;</FONT></B>) )
      online_mode = command_line.next(online_mode);
  
    Mesh mesh(init.comm());
  
    <B><FONT COLOR="#5F9EA0">MeshTools</FONT></B>::Generation::build_cube (mesh,
                                       n_elem_xy, n_elem_xy, n_elem_z,
                                       -0.2, 0.2,
                                       -0.2, 0.2,
                                       0., 3.,
                                       HEX8);
  
    EquationSystems equation_systems (mesh);
  
    SimpleEIMConstruction &amp; eim_construction =
      equation_systems.add_system&lt;SimpleEIMConstruction&gt; (<B><FONT COLOR="#BC8F8F">&quot;EIM&quot;</FONT></B>);
    SimpleRBConstruction &amp; rb_construction =
      equation_systems.add_system&lt;SimpleRBConstruction&gt; (<B><FONT COLOR="#BC8F8F">&quot;RB&quot;</FONT></B>);
  
    equation_systems.init ();
  
    equation_systems.print_info();
    mesh.print_info();
  
    SimpleRBEvaluation rb_eval(mesh.comm());
  
    SimpleEIMEvaluation eim_rb_eval(mesh.comm());
  
    eim_construction.set_rb_evaluation(eim_rb_eval);
    rb_construction.set_rb_evaluation(rb_eval);
  
    <B><FONT COLOR="#A020F0">if</FONT></B>(!online_mode) <I><FONT COLOR="#B22222">// Perform the Offline stage of the RB method
</FONT></I>    {
      eim_construction.process_parameters_file(eim_parameters);
      eim_construction.print_info();
  
      eim_construction.initialize_rb_construction();
  
      eim_construction.train_reduced_basis();
      eim_construction.get_rb_evaluation().write_offline_data_to_files(<B><FONT COLOR="#BC8F8F">&quot;eim_data&quot;</FONT></B>);
  
      rb_construction.process_parameters_file(rb_parameters);
  
      eim_rb_eval.initialize_eim_theta_objects();
      rb_eval.get_rb_theta_expansion().attach_multiple_A_theta(eim_rb_eval.get_eim_theta_objects());
  
      eim_construction.initialize_eim_assembly_objects();
      rb_construction.get_rb_assembly_expansion().attach_multiple_A_assembly(eim_construction.get_eim_assembly_objects());
  
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
    <B><FONT COLOR="#A020F0">else</FONT></B> <I><FONT COLOR="#B22222">// Perform the Online stage of the RB method
</FONT></I>    {
      eim_rb_eval.read_offline_data_from_files(<B><FONT COLOR="#BC8F8F">&quot;eim_data&quot;</FONT></B>);
  
      eim_rb_eval.initialize_eim_theta_objects();
      rb_eval.get_rb_theta_expansion().attach_multiple_A_theta(eim_rb_eval.get_eim_theta_objects());
  
      rb_eval.read_offline_data_from_files(<B><FONT COLOR="#BC8F8F">&quot;rb_data&quot;</FONT></B>);
  
      Real online_curvature = infile(<B><FONT COLOR="#BC8F8F">&quot;online_curvature&quot;</FONT></B>, 0.);
      Real online_Bi        = infile(<B><FONT COLOR="#BC8F8F">&quot;online_Bi&quot;</FONT></B>, 0.);
      Real online_kappa     = infile(<B><FONT COLOR="#BC8F8F">&quot;online_kappa&quot;</FONT></B>, 0.);
      RBParameters online_mu;
      online_mu.set_value(<B><FONT COLOR="#BC8F8F">&quot;curvature&quot;</FONT></B>, online_curvature);
      online_mu.set_value(<B><FONT COLOR="#BC8F8F">&quot;Bi&quot;</FONT></B>, online_Bi);
      online_mu.set_value(<B><FONT COLOR="#BC8F8F">&quot;kappa&quot;</FONT></B>, online_kappa);
      rb_eval.set_parameters(online_mu);
      rb_eval.print_parameters();
      rb_eval.rb_solve( rb_eval.get_n_basis_functions() );
  
      <B><FONT COLOR="#A020F0">if</FONT></B>(store_basis_functions)
      {
        eim_rb_eval.read_in_basis_functions(eim_construction,<B><FONT COLOR="#BC8F8F">&quot;eim_data&quot;</FONT></B>);
        rb_eval.read_in_basis_functions(rb_construction,<B><FONT COLOR="#BC8F8F">&quot;rb_data&quot;</FONT></B>);
  
        eim_construction.load_rb_solution();
        rb_construction.load_rb_solution();
  
        transform_mesh_and_plot(equation_systems,online_curvature,<B><FONT COLOR="#BC8F8F">&quot;RB_sol.e&quot;</FONT></B>);
      }
    }
  
    <B><FONT COLOR="#A020F0">return</FONT></B> 0;
  }
  
  <B><FONT COLOR="#228B22">void</FONT></B> transform_mesh_and_plot(EquationSystems&amp; es, Real curvature, <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp; filename)
  {
    MeshBase&amp; mesh = es.get_mesh();
  
    <B><FONT COLOR="#5F9EA0">MeshBase</FONT></B>::node_iterator       node_it  = mesh.nodes_begin();
    <B><FONT COLOR="#228B22">const</FONT></B> MeshBase::node_iterator node_end = mesh.nodes_end();
  
    <B><FONT COLOR="#A020F0">for</FONT></B>( ; node_it != node_end; node_it++)
    {
      Node* node = *node_it;
  
      Real x = (*node)(0);
      Real z = (*node)(2);
  
      (*node)(0) = -1./curvature + (1./curvature + x)*cos(curvature*z);
      (*node)(2) = (1./curvature + x)*sin(curvature*z);
    }
  
  #ifdef LIBMESH_HAVE_EXODUS_API
    ExodusII_IO(mesh).write_equation_systems(filename, es);
  #endif
  }
</pre> 
<a name="output"></a> 
<br><br><br> <h1> The console output of the program: </h1> 
<pre>
make[4]: Entering directory `/net/spark/workspace/roystgnr/libmesh/git/devel/examples/reduced_basis/reduced_basis_ex6'
***************************************************************
* Running Example reduced_basis_ex6:
*  mpirun -np 4 example-devel -online_mode 0 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc
***************************************************************
 
*** Warning, This code is untested, experimental, or likely to see future API changes: ../src/reduced_basis/rb_parametrized.C, line 41, compiled Apr 19 2013 at 11:42:51 ***
 EquationSystems
  n_systems()=2
   System #0, "EIM"
    Type "RBConstruction"
    Variables={ "x_comp_of_G" "y_comp_of_G" "z_comp_of_G" } 
    Finite Element Types="LAGRANGE", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" 
    Approximation Orders="FIRST", "THIRD" 
    n_dofs()=18513
    n_local_dofs()=4926
    n_constrained_dofs()=0
    n_local_constrained_dofs()=0
    n_vectors()=1
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 23.06
      Average Off-Processor Bandwidth <= 0.964511
      Maximum  On-Processor Bandwidth <= 33
      Maximum Off-Processor Bandwidth <= 14
    DofMap Constraints
      Number of DoF Constraints = 0
      Number of Node Constraints = 0
   System #1, "RB"
    Type "RBConstruction"
    Variables="u" 
    Finite Element Types="LAGRANGE", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" 
    Approximation Orders="FIRST", "THIRD" 
    n_dofs()=6171
    n_local_dofs()=1642
    n_constrained_dofs()=242
    n_local_constrained_dofs()=121
    n_vectors()=1
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 23.06
      Average Off-Processor Bandwidth <= 0.964511
      Maximum  On-Processor Bandwidth <= 33
      Maximum Off-Processor Bandwidth <= 14
    DofMap Constraints
      Number of DoF Constraints = 242
      Average DoF Constraint Length= 0
      Number of Node Constraints = 0

 Mesh Information:
  mesh_dimension()=3
  spatial_dimension()=3
  n_nodes()=6171
    n_local_nodes()=1642
  n_elem()=5000
    n_local_elem()=1250
    n_active_elem()=5000
  n_subdomains()=1
  n_partitions()=4
  n_processors()=4
  n_threads()=1
  processor_id()=0

Initializing training parameters with deterministic training set...
Parameter curvature: log scaling = 0


RBConstruction parameters:
system name: EIM
constrained_problem: 0
Nmax: 20
Basis training error tolerance: 0.001
Aq operators attached: 0
Fq functions attached: 0
n_outputs: 0
Number of parameters: 1
Parameter curvature: Min = 0.1, Max = 1.0472, value = 0.5
n_training_samples: 25
single-matrix mode? 0
reuse preconditioner? 1
use a relative error bound in greedy? 1
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
curvature: 1.047200e+00

Enriching the RB space
Updating RB matrices

---- Basis dimension: 1 ----
Performing RB solves on training set
Maximum (relative) error bound is 0.221031

Performing truth solve at parameter:
curvature: 1.000000e-01

Enriching the RB space
Updating RB matrices

---- Basis dimension: 2 ----
Performing RB solves on training set
Maximum (relative) error bound is 0.0106172

Performing truth solve at parameter:
curvature: 6.130667e-01

Enriching the RB space
Updating RB matrices

---- Basis dimension: 3 ----
Performing RB solves on training set
Maximum (relative) error bound is 0.000307645

Specified error tolerance reached.
Perform one more Greedy iteration for error bounds.
Performing truth solve at parameter:
curvature: 2.973333e-01

Enriching the RB space
Updating RB matrices

---- Basis dimension: 3 ----
Performing RB solves on training set
Maximum (relative) error bound is 0.000307645

Extra Greedy iteration finished.
Initializing training parameters with random training set...
Parameter Bi: log scaling = 0
Parameter curvature: log scaling = 1
Parameter kappa: log scaling = 1


RBConstruction parameters:
system name: RB
constrained_problem: 0
Nmax: 15
Basis training error tolerance: 0.001
Aq operators attached: 6
Fq functions attached: 2
n_outputs: 0
Number of parameters: 3
Parameter Bi: Min = 0.001, Max = 0.01, value = 0.005
Parameter curvature: Min = 0.1, Max = 1.0472, value = 0.5
Parameter kappa: Min = 0.5, Max = 2, value = 1
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
Bi: 2.202196e-03
curvature: 1.868130e-01
kappa: 9.085363e-01

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 1 ----
Performing RB solves on training set
Maximum (relative) error bound is 0.0529974

Performing truth solve at parameter:
Bi: 2.188696e-03
curvature: 1.028177e+00
kappa: 1.916434e+00

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 2 ----
Performing RB solves on training set
Maximum (relative) error bound is 0.0182385

Performing truth solve at parameter:
Bi: 9.883868e-03
curvature: 2.730655e-01
kappa: 1.928031e+00

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 3 ----
Performing RB solves on training set
Maximum (relative) error bound is 0.00101092

Performing truth solve at parameter:
Bi: 9.493733e-03
curvature: 9.310653e-01
kappa: 1.824189e+00

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 4 ----
Performing RB solves on training set
Maximum (relative) error bound is 0.000949469

Specified error tolerance reached.

 -------------------------------------------------------------------------------------------------------------------
| Processor id:   0                                                                                                 |
| Num Processors: 4                                                                                                 |
| Time:           Fri Apr 19 11:57:21 2013                                                                          |
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
| libMesh Performance: Alive time=15.6312, Active time=15.4082                                                      |
 -------------------------------------------------------------------------------------------------------------------
| Event                                 nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                                 w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|-------------------------------------------------------------------------------------------------------------------|
|                                                                                                                   |
|                                                                                                                   |
| DofMap                                                                                                            |
|   add_neighbors_to_send_list()        2         0.0184      0.009177    0.0224      0.011218    0.12     0.15     |
|   build_constraint_matrix()           11250     0.0178      0.000002    0.0178      0.000002    0.12     0.12     |
|   build_sparsity()                    2         0.0178      0.008877    0.0644      0.032211    0.12     0.42     |
|   cnstrn_elem_mat_vec()               11250     0.0207      0.000002    0.0207      0.000002    0.13     0.13     |
|   create_dof_constraints()            2         0.0125      0.006243    0.0372      0.018624    0.08     0.24     |
|   distribute_dofs()                   2         0.0159      0.007932    0.0805      0.040237    0.10     0.52     |
|   dof_indices()                       194228    1.9526      0.000010    1.9526      0.000010    12.67    12.67    |
|   prepare_send_list()                 2         0.0001      0.000075    0.0001      0.000075    0.00     0.00     |
|   reinit()                            2         0.0301      0.015034    0.0301      0.015034    0.20     0.20     |
|                                                                                                                   |
| FE                                                                                                                |
|   compute_shape_functions()           109420    1.7810      0.000016    1.7810      0.000016    11.56    11.56    |
|   init_shape_functions()              11998     0.1760      0.000015    0.1760      0.000015    1.14     1.14     |
|   inverse_map()                       90030     0.6992      0.000008    0.6992      0.000008    4.54     4.54     |
|                                                                                                                   |
| FEMap                                                                                                             |
|   compute_affine_map()                109420    0.6651      0.000006    0.6651      0.000006    4.32     4.32     |
|   compute_face_map()                  11920     0.0519      0.000004    0.0519      0.000004    0.34     0.34     |
|   init_face_shape_functions()         20        0.0002      0.000010    0.0002      0.000010    0.00     0.00     |
|   init_reference_to_physical_map()    11998     0.1490      0.000012    0.1490      0.000012    0.97     0.97     |
|                                                                                                                   |
| Mesh                                                                                                              |
|   find_neighbors()                    1         0.0133      0.013326    0.0288      0.028794    0.09     0.19     |
|   renumber_nodes_and_elem()           2         0.0016      0.000775    0.0016      0.000775    0.01     0.01     |
|                                                                                                                   |
| MeshCommunication                                                                                                 |
|   assign_global_indices()             2         0.2356      0.117802    0.2376      0.118812    1.53     1.54     |
|   compute_hilbert_indices()           2         0.0268      0.013388    0.0268      0.013388    0.17     0.17     |
|   find_global_indices()               2         0.0038      0.001888    0.0586      0.029282    0.02     0.38     |
|   parallel_sort()                     2         0.0009      0.000438    0.0264      0.013194    0.01     0.17     |
|                                                                                                                   |
| MeshTools::Generation                                                                                             |
|   build_cube()                        1         0.0034      0.003354    0.0034      0.003354    0.02     0.02     |
|                                                                                                                   |
| MetisPartitioner                                                                                                  |
|   partition()                         1         0.0972      0.097236    0.1290      0.129030    0.63     0.84     |
|                                                                                                                   |
| Parallel                                                                                                          |
|   allgather()                         22        0.0277      0.001261    0.0281      0.001278    0.18     0.18     |
|   barrier()                           2         0.0000      0.000011    0.0000      0.000011    0.00     0.00     |
|   broadcast()                         45        0.0002      0.000005    0.0002      0.000005    0.00     0.00     |
|   max(bool)                           10        0.0001      0.000006    0.0001      0.000006    0.00     0.00     |
|   max(scalar)                         206       0.0023      0.000011    0.0023      0.000011    0.02     0.02     |
|   max(vector)                         45        0.0003      0.000007    0.0008      0.000018    0.00     0.01     |
|   maxloc(scalar)                      13        0.1915      0.014734    0.1915      0.014734    1.24     1.24     |
|   min(bool)                           220       0.0008      0.000004    0.0008      0.000004    0.01     0.01     |
|   min(scalar)                         180       0.1381      0.000767    0.1381      0.000767    0.90     0.90     |
|   min(vector)                         45        0.0004      0.000009    0.0039      0.000087    0.00     0.03     |
|   probe()                             142       0.0070      0.000049    0.0070      0.000049    0.05     0.05     |
|   receive()                           110       0.0011      0.000010    0.0079      0.000071    0.01     0.05     |
|   send()                              86        0.0003      0.000004    0.0003      0.000004    0.00     0.00     |
|   send_receive()                      90        0.0005      0.000006    0.0081      0.000090    0.00     0.05     |
|   sum()                               45        0.0257      0.000572    0.0258      0.000574    0.17     0.17     |
|                                                                                                                   |
| Parallel::Request                                                                                                 |
|   wait()                              86        0.0001      0.000001    0.0001      0.000001    0.00     0.00     |
|                                                                                                                   |
| Partitioner                                                                                                       |
|   set_node_processor_ids()            1         0.0026      0.002606    0.0471      0.047085    0.02     0.31     |
|   set_parent_processor_ids()          1         0.0007      0.000725    0.0007      0.000725    0.00     0.00     |
|                                                                                                                   |
| PetscLinearSolver                                                                                                 |
|   solve()                             55        4.2689      0.077616    4.2689      0.077616    27.71    27.71    |
|                                                                                                                   |
| PointLocatorTree                                                                                                  |
|   init(no master)                     1         0.0139      0.013881    0.0146      0.014632    0.09     0.09     |
|   operator()                          15        0.0003      0.000021    0.0005      0.000035    0.00     0.00     |
|                                                                                                                   |
| RBConstruction                                                                                                    |
|   add_scaled_matrix_and_vector()      10        1.2825      0.128252    3.4498      0.344981    8.32     22.39    |
|   clear()                             3         0.0013      0.000421    0.0013      0.000421    0.01     0.01     |
|   compute_Fq_representor_innerprods() 2         0.0170      0.008522    0.2126      0.106319    0.11     1.38     |
|   compute_max_error_bound()           9         0.0071      0.000785    0.5623      0.062476    0.05     3.65     |
|   enrich_RB_space()                   4         0.0024      0.000596    0.0024      0.000596    0.02     0.02     |
|   train_reduced_basis()               2         0.0022      0.001098    3.7907      1.895364    0.01     24.60    |
|   truth_assembly()                    4         0.0954      0.023853    0.0955      0.023886    0.62     0.62     |
|   truth_solve()                       4         0.0013      0.000328    0.6507      0.162687    0.01     4.22     |
|   update_RB_system_matrices()         8         0.0268      0.003347    0.0268      0.003347    0.17     0.17     |
|   update_residual_terms()             4         0.0871      0.021780    1.7305      0.432620    0.57     11.23    |
|                                                                                                                   |
| RBEIMConstruction                                                                                                 |
|   compute_best_fit_error()            100       0.1345      0.001345    0.1446      0.001446    0.87     0.94     |
|   enrich_RB_space()                   4         0.1093      0.027313    0.6010      0.150245    0.71     3.90     |
|   truth_solve()                       129       2.6739      0.020728    7.4348      0.057634    17.35    48.25    |
|   update_RB_system_matrices()         4         0.0008      0.000192    0.0075      0.001874    0.00     0.05     |
|                                                                                                                   |
| RBEIMEvaluation                                                                                                   |
|   rb_solve()                          1256      0.0055      0.000004    0.0055      0.000004    0.04     0.04     |
|   write_offline_data_to_files()       1         0.0004      0.000362    0.0011      0.001070    0.00     0.01     |
|                                                                                                                   |
| RBEvaluation                                                                                                      |
|   clear()                             3         0.0001      0.000037    0.0001      0.000037    0.00     0.00     |
|   compute_residual_dual_norm()        1250      0.2195      0.000176    0.2195      0.000176    1.42     1.42     |
|   rb_solve()                          1250      0.0129      0.000010    0.2380      0.000190    0.08     1.54     |
|   resize_data_structures()            2         0.0001      0.000041    0.0001      0.000041    0.00     0.00     |
|   write_offline_data_to_files()       2         0.0019      0.000933    0.0019      0.000933    0.01     0.01     |
|   write_out_basis_functions()         2         0.0006      0.000276    0.2931      0.146554    0.00     1.90     |
|   write_out_vectors()                 2         0.0524      0.026214    0.2926      0.146278    0.34     1.90     |
 -------------------------------------------------------------------------------------------------------------------
| Totals:                               567034    15.4082                                         100.00            |
 -------------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example reduced_basis_ex6:
*  mpirun -np 4 example-devel -online_mode 0 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc
***************************************************************
***************************************************************
* Running Example reduced_basis_ex6:
*  mpirun -np 4 example-devel -online_mode 1 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc
***************************************************************
 
*** Warning, This code is untested, experimental, or likely to see future API changes: ../src/reduced_basis/rb_parametrized.C, line 41, compiled Apr 19 2013 at 11:42:51 ***
 EquationSystems
  n_systems()=2
   System #0, "EIM"
    Type "RBConstruction"
    Variables={ "x_comp_of_G" "y_comp_of_G" "z_comp_of_G" } 
    Finite Element Types="LAGRANGE", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" 
    Approximation Orders="FIRST", "THIRD" 
    n_dofs()=18513
    n_local_dofs()=4926
    n_constrained_dofs()=0
    n_local_constrained_dofs()=0
    n_vectors()=1
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 23.06
      Average Off-Processor Bandwidth <= 0.964511
      Maximum  On-Processor Bandwidth <= 33
      Maximum Off-Processor Bandwidth <= 14
    DofMap Constraints
      Number of DoF Constraints = 0
      Number of Node Constraints = 0
   System #1, "RB"
    Type "RBConstruction"
    Variables="u" 
    Finite Element Types="LAGRANGE", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" 
    Approximation Orders="FIRST", "THIRD" 
    n_dofs()=6171
    n_local_dofs()=1642
    n_constrained_dofs()=242
    n_local_constrained_dofs()=121
    n_vectors()=1
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 23.06
      Average Off-Processor Bandwidth <= 0.964511
      Maximum  On-Processor Bandwidth <= 33
      Maximum Off-Processor Bandwidth <= 14
    DofMap Constraints
      Number of DoF Constraints = 242
      Average DoF Constraint Length= 0
      Number of Node Constraints = 0

 Mesh Information:
  mesh_dimension()=3
  spatial_dimension()=3
  n_nodes()=6171
    n_local_nodes()=1642
  n_elem()=5000
    n_local_elem()=1250
    n_active_elem()=5000
  n_subdomains()=1
  n_partitions()=4
  n_processors()=4
  n_threads()=1
  processor_id()=0

Bi: 5.000000e-03
curvature: 1.047200e+00
kappa: 1.300000e+00


 -------------------------------------------------------------------------------------------------------------------
| Processor id:   0                                                                                                 |
| Num Processors: 4                                                                                                 |
| Time:           Fri Apr 19 11:57:22 2013                                                                          |
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
| libMesh Performance: Alive time=1.03509, Active time=0.999776                                                |
 --------------------------------------------------------------------------------------------------------------
| Event                            nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                            w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|--------------------------------------------------------------------------------------------------------------|
|                                                                                                              |
|                                                                                                              |
| DofMap                                                                                                       |
|   add_neighbors_to_send_list()   2         0.0344      0.017207    0.0419      0.020955    3.44     4.19     |
|   build_sparsity()               2         0.0321      0.016037    0.0764      0.038178    3.21     7.64     |
|   create_dof_constraints()       2         0.0209      0.010470    0.0663      0.033174    2.09     6.64     |
|   distribute_dofs()              2         0.0297      0.014864    0.0885      0.044226    2.97     8.85     |
|   dof_indices()                  15478     0.1173      0.000008    0.1173      0.000008    11.73    11.73    |
|   prepare_send_list()            2         0.0002      0.000088    0.0002      0.000088    0.02     0.02     |
|   reinit()                       2         0.0579      0.028949    0.0579      0.028949    5.79     5.79     |
|                                                                                                              |
| EquationSystems                                                                                              |
|   build_solution_vector()        1         0.0069      0.006854    0.0375      0.037517    0.69     3.75     |
|                                                                                                              |
| ExodusII_IO                                                                                                  |
|   write_nodal_data()             1         0.1643      0.164343    0.1643      0.164343    16.44    16.44    |
|                                                                                                              |
| Mesh                                                                                                         |
|   find_neighbors()               1         0.0265      0.026501    0.0266      0.026596    2.65     2.66     |
|   renumber_nodes_and_elem()      2         0.0024      0.001223    0.0024      0.001223    0.24     0.24     |
|                                                                                                              |
| MeshCommunication                                                                                            |
|   assign_global_indices()        2         0.1856      0.092795    0.1986      0.099311    18.56    19.87    |
|   compute_hilbert_indices()      2         0.0460      0.022997    0.0460      0.022997    4.60     4.60     |
|   find_global_indices()          2         0.0062      0.003085    0.0554      0.027689    0.62     5.54     |
|   parallel_sort()                2         0.0011      0.000529    0.0028      0.001402    0.11     0.28     |
|                                                                                                              |
| MeshOutput                                                                                                   |
|   write_equation_systems()       1         0.0001      0.000071    0.2020      0.202007    0.01     20.21    |
|                                                                                                              |
| MeshTools::Generation                                                                                        |
|   build_cube()                   1         0.0058      0.005840    0.0058      0.005840    0.58     0.58     |
|                                                                                                              |
| MetisPartitioner                                                                                             |
|   partition()                    1         0.1852      0.185179    0.2128      0.212833    18.52    21.29    |
|                                                                                                              |
| Parallel                                                                                                     |
|   allgather()                    22        0.0111      0.000505    0.0112      0.000508    1.11     1.12     |
|   barrier()                      2         0.0000      0.000015    0.0000      0.000015    0.00     0.00     |
|   broadcast()                    48        0.0001      0.000003    0.0001      0.000002    0.01     0.01     |
|   max(bool)                      2         0.0000      0.000004    0.0000      0.000004    0.00     0.00     |
|   max(scalar)                    177       0.0009      0.000005    0.0009      0.000005    0.09     0.09     |
|   max(vector)                    39        0.0003      0.000007    0.0007      0.000018    0.03     0.07     |
|   min(bool)                      201       0.0008      0.000004    0.0008      0.000004    0.08     0.08     |
|   min(scalar)                    166       0.0144      0.000087    0.0144      0.000087    1.44     1.44     |
|   min(vector)                    39        0.0003      0.000009    0.0011      0.000029    0.03     0.11     |
|   probe()                        110       0.0012      0.000010    0.0012      0.000010    0.12     0.12     |
|   receive()                      98        0.0006      0.000006    0.0017      0.000017    0.06     0.17     |
|   send()                         98        0.0004      0.000004    0.0004      0.000004    0.04     0.04     |
|   send_receive()                 90        0.0006      0.000007    0.0026      0.000029    0.06     0.26     |
|   sum()                          48        0.0040      0.000084    0.0131      0.000273    0.40     1.31     |
|                                                                                                              |
| Parallel::Request                                                                                            |
|   wait()                         98        0.0004      0.000005    0.0004      0.000005    0.04     0.04     |
|                                                                                                              |
| Partitioner                                                                                                  |
|   set_node_processor_ids()       1         0.0047      0.004702    0.0049      0.004870    0.47     0.49     |
|   set_parent_processor_ids()     1         0.0016      0.001567    0.0016      0.001567    0.16     0.16     |
|                                                                                                              |
| RBConstruction                                                                                               |
|   clear()                        3         0.0004      0.000137    0.0004      0.000137    0.04     0.04     |
|   load_rb_solution()             2         0.0003      0.000140    0.0003      0.000140    0.03     0.03     |
|                                                                                                              |
| RBEIMEvaluation                                                                                              |
|   rb_solve()                     1         0.0064      0.006352    0.0064      0.006352    0.64     0.64     |
|   read_offline_data_from_files() 1         0.0001      0.000124    0.0005      0.000453    0.01     0.05     |
|                                                                                                              |
| RBEvaluation                                                                                                 |
|   clear()                        3         0.0001      0.000019    0.0001      0.000019    0.01     0.01     |
|   compute_residual_dual_norm()   1         0.0007      0.000727    0.0007      0.000727    0.07     0.07     |
|   rb_solve()                     1         0.0001      0.000081    0.0072      0.007160    0.01     0.72     |
|   read_in_basis_functions()      2         0.0000      0.000017    0.2334      0.116703    0.00     23.35    |
|   read_in_vectors()              2         0.0270      0.013476    0.2334      0.116686    2.70     23.34    |
|   read_offline_data_from_files() 2         0.0006      0.000310    0.0007      0.000339    0.06     0.07     |
|   resize_data_structures()       2         0.0001      0.000028    0.0001      0.000028    0.01     0.01     |
 --------------------------------------------------------------------------------------------------------------
| Totals:                          16766     0.9998                                          100.00            |
 --------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example reduced_basis_ex6:
*  mpirun -np 4 example-devel -online_mode 1 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc
***************************************************************
make[4]: Leaving directory `/net/spark/workspace/roystgnr/libmesh/git/devel/examples/reduced_basis/reduced_basis_ex6'
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
