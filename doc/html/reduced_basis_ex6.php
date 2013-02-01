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
        
          SimpleEIMEvaluation()
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
        
        #endif</pre>
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
                                const std::string& name,
                                const unsigned int number)
          : Parent(es, name, number),
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
<h1>Reduced Basis: Example 6 - Heat transfer on a curved domain in 3D</h1>


<br><br>In this example we consider heat transfer modeled by a Poisson equation with
Robin boundary condition:
-kappa \Laplacian u = 1, on \Omega
-kappa du\dn = kappa Bi u, on \partial\Omega_Biot,
u = 0 on \partial\Omega_Dirichlet,

<br><br>We consider a reference domain \Omega_hat = [-0.2,0.2]x[-0.2,0.2]x[0,3], and the
physical domain is then obtain via the parametrized mapping:
x = -1/mu + (1/mu+x_hat)*cos(mu*z_hat)
y = y_hat
z = (1/mu+x_hat)*sin(mu*z_hat)
for (x_hat,y_hat,z_hat) \in \Omega_hat. (Here "hats" denotes reference domain.)
Also, the "reference Dirichlet boundaries" are [-0.2,0.2]x[-0.2,0.2]x{0} and
[-0.2,0.2]x[-0.2,0.2]x{3}, and the remaining boundaries are the "Biot" boundaries.


<br><br>Then, after putting the PDE into weak form and mapping it to the reference domain,
we obtain:
\kappa \int_\Omega_hat [ (1+mu*x_hat) v_x w_x + (1+mu*x_hat) v_y w_y + 1/(1+mu*x_hat) v_z w_z ]
+ \kappa Bi \int_\partial\Omega_hat_Biot1 (1-0.2mu) u v
+ \kappa Bi \int_\partial\Omega_hat_Biot2 (1+mu x_hat) u v
+ \kappa Bi \int_\partial\Omega_hat_Biot3 (1+0.2mu) u v
= \int_\Omega_hat (1+mu x_hat) v
where
\partial\Omega_hat_Biot1 = [-0.2] x [-0.2,0.2] x [0,3]
\partial\Omega_hat_Biot2 = [-0.2,0.2] x {-0.2} x [0,3] \UNION [-0.2,0.2] x {0.2} x [0,3]
\partial\Omega_hat_Biot3 = [0.2] x [-0.2,0.2] x [0,3]


<br><br>The term
\kappa \int_\Omega_hat 1/(1+mu*x_hat) v_z w_z 
is "non-affine" (in the Reduced Basis sense), since we can't express it
in the form \sum theta_q(kappa,mu) a(v,w). As a result, (as in
reduced_basis_ex4) we must employ the Empirical Interpolation Method (EIM)
in order to apply the Reduced Basis method here.


<br><br>The approach we use is to construct an EIM approximation, G_EIM, to the vector-valued function
G(x_hat,y_hat;mu) = (1 + mu*x_hat, 1 + mu*x_hat, 1/(1+mu*x_hat))
and then we express the "volumetric integral part" of the left-hand side operator as
a(v,w;mu) = \int_\hat\Omega G_EIM(x_hat,y_hat;mu) \dot (v_x w_x, v_y w_y, v_z w_z).
(We actually only need EIM for the third component of G_EIM, but it's helpful to
demonstrate "vector-valued" EIM here.)


<br><br>

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
Build a mesh.
</div>

<div class ="fragment">
<pre>
          Mesh mesh;
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
  
    SimpleEIMEvaluation()
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
                          <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp; name,
                          <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> number)
    : Parent(es, name, number),
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
  
    Mesh mesh;
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
  
    SimpleRBEvaluation rb_eval;
  
    SimpleEIMEvaluation eim_rb_eval;
    
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
***************************************************************
* Running Example reduced_basis_ex6:
*  mpirun -np 12 example-devel -online_mode 0 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
***************************************************************
 
*** Warning, This code is untested, experimental, or likely to see future API changes: src/reduced_basis/rb_parametrized.C, line 41, compiled Jan 31 2013 at 21:51:32 ***
 EquationSystems
  n_systems()=2
   System #0, "EIM"
    Type "RBConstruction"
    Variables={ "x_comp_of_G" "y_comp_of_G" "z_comp_of_G" } 
    Finite Element Types="LAGRANGE", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" 
    Approximation Orders="FIRST", "THIRD" 
    n_dofs()=18513
    n_local_dofs()=1911
    n_constrained_dofs()=0
    n_local_constrained_dofs()=0
    n_vectors()=1
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 22.0431
      Average Off-Processor Bandwidth <= 3.38324
      Maximum  On-Processor Bandwidth <= 39
      Maximum Off-Processor Bandwidth <= 27
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
    n_local_dofs()=637
    n_constrained_dofs()=242
    n_local_constrained_dofs()=0
    n_vectors()=1
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 22.0431
      Average Off-Processor Bandwidth <= 3.38324
      Maximum  On-Processor Bandwidth <= 39
      Maximum Off-Processor Bandwidth <= 27
    DofMap Constraints
      Number of DoF Constraints = 242
      Average DoF Constraint Length= 0
      Number of Node Constraints = 0

 Mesh Information:
  mesh_dimension()=3
  spatial_dimension()=3
  n_nodes()=6171
    n_local_nodes()=637
  n_elem()=5000
    n_local_elem()=424
    n_active_elem()=5000
  n_subdomains()=1
  n_partitions()=12
  n_processors()=12
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
Bi: 4.520254e-03
curvature: 1.159800e-01
kappa: 5.334842e-01

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 1 ----
Performing RB solves on training set
Maximum (relative) error bound is 0.0553459

Performing truth solve at parameter:
Bi: 9.542854e-03
curvature: 9.665050e-01
kappa: 1.966371e+00

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 2 ----
Performing RB solves on training set
Maximum (relative) error bound is 0.0169562

Performing truth solve at parameter:
Bi: 2.093504e-03
curvature: 8.789475e-01
kappa: 1.977278e+00

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 3 ----
Performing RB solves on training set
Maximum (relative) error bound is 0.00120136

Performing truth solve at parameter:
Bi: 9.664390e-03
curvature: 5.546889e-01
kappa: 1.948182e+00

Enriching the RB space
Updating RB matrices
Updating RB residual terms

---- Basis dimension: 4 ----
Performing RB solves on training set
Maximum (relative) error bound is 0.00077215

Specified error tolerance reached.
************************************************************************************************************************
***             WIDEN YOUR WINDOW TO 120 CHARACTERS.  Use 'enscript -r -fCourier9' to print this document            ***
************************************************************************************************************************

---------------------------------------------- PETSc Performance Summary: ----------------------------------------------

/workspace/libmesh/examples/reduced_basis/reduced_basis_ex6/.libs/lt-example-devel on a intel-12. named hbar.ices.utexas.edu with 12 processors, by benkirk Thu Jan 31 22:20:06 2013
Using Petsc Release Version 3.3.0, Patch 2, Fri Jul 13 15:42:00 CDT 2012 

                         Max       Max/Min        Avg      Total 
Time (sec):           2.594e+01      1.00005   2.594e+01
Objects:              6.260e+02      1.00000   6.260e+02
Flops:                1.521e+09      2.51341   9.747e+08  1.170e+10
Flops/sec:            5.865e+07      2.51341   3.758e+07  4.509e+08
MPI Messages:         4.266e+04      3.50000   2.945e+04  3.534e+05
MPI Message Lengths:  2.237e+07      1.98892   6.093e+02  2.154e+08
MPI Reductions:       1.477e+04      1.00000

Flop counting convention: 1 flop = 1 real number operation of type (multiply/divide/add/subtract)
                            e.g., VecAXPY() for real vectors of length N --> 2N flops
                            and VecAXPY() for complex vectors of length N --> 8N flops

Summary of Stages:   ----- Time ------  ----- Flops -----  --- Messages ---  -- Message Lengths --  -- Reductions --
                        Avg     %Total     Avg     %Total   counts   %Total     Avg         %Total   counts   %Total 
 0:      Main Stage: 2.5940e+01 100.0%  1.1697e+10 100.0%  3.534e+05 100.0%  6.093e+02      100.0%  1.477e+04 100.0% 

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

VecDot               739 1.0 2.9132e-02 2.3 1.58e+06 1.7 0.0e+00 0.0e+00 7.4e+02  0  0  0  0  5   0  0  0  0  5   511
VecMDot             4889 1.0 7.1474e-01 5.6 1.08e+08 1.7 0.0e+00 0.0e+00 4.9e+03  2  9  0  0 33   2  9  0  0 33  1418
VecNorm             5228 1.0 2.5650e-01 4.0 8.17e+06 1.7 0.0e+00 0.0e+00 5.2e+03  1  1  0  0 35   1  1  0  0 35   300
VecScale            5113 1.0 6.2265e-03 1.7 4.05e+06 1.7 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0  6122
VecCopy              539 1.0 1.4906e-03 1.4 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecSet              6327 1.0 8.2526e-03 1.3 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecAXPY              615 1.0 1.2749e-02 1.1 1.54e+06 1.7 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0  1137
VecMAXPY            5073 1.0 6.8485e-02 1.7 1.15e+08 1.7 0.0e+00 0.0e+00 0.0e+00  0  9  0  0  0   0  9  0  0  0 15880
VecAssemblyBegin     475 1.0 1.8884e-0112.7 0.00e+00 0.0 1.5e+03 3.7e+03 1.4e+03  1  0  0  3 10   1  0  0  3 10     0
VecAssemblyEnd       475 1.0 7.5340e-04 2.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecScatterBegin     5966 1.0 4.2655e-02 1.7 0.00e+00 0.0 3.5e+05 5.9e+02 0.0e+00  0  0 98 94  0   0  0 98 94  0     0
VecScatterEnd       5966 1.0 1.9933e-0114.3 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  1  0  0  0  0   1  0  0  0  0     0
VecNormalize        5073 1.0 2.5043e-01 3.7 1.21e+07 1.7 0.0e+00 0.0e+00 5.1e+03  1  1  0  0 34   1  1  0  0 34   453
MatMult             5073 1.0 4.0862e-01 1.6 1.84e+08 1.7 2.9e+05 5.2e+02 0.0e+00  1 15 83 72  0   1 15 83 72  0  4263
MatMultAdd           714 1.0 4.4505e-02 1.4 3.63e+07 1.7 4.1e+04 7.2e+02 0.0e+00  0  3 12 14  0   0  3 12 14  0  7708
MatSolve            5128 1.0 9.9276e-01 2.6 1.01e+09 3.1 0.0e+00 0.0e+00 0.0e+00  2 60  0  0  0   2 60  0  0  0  7077
MatLUFactorNum        10 1.0 7.0593e-02 4.7 5.93e+07 5.3 0.0e+00 0.0e+00 0.0e+00  0  3  0  0  0   0  3  0  0  0  4824
MatILUFactorSym       10 1.0 1.9350e-01 4.7 0.00e+00 0.0 0.0e+00 0.0e+00 3.0e+01  0  0  0  0  0   0  0  0  0  0     0
MatAssemblyBegin     874 1.0 1.5920e+0047.9 0.00e+00 0.0 6.7e+02 9.1e+03 1.7e+03  4  0  0  3 12   4  0  0  3 12     0
MatAssemblyEnd       874 1.0 4.7589e-02 1.3 0.00e+00 0.0 4.5e+03 1.3e+02 3.3e+02  0  0  1  0  2   0  0  1  0  2     0
MatGetRow          45850 1.7 9.2087e-03 1.5 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatGetRowIJ           10 1.0 1.5736e-05 5.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatGetOrdering        10 1.0 4.0722e-04 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 2.0e+01  0  0  0  0  0   0  0  0  0  0     0
MatZeroEntries        39 1.0 1.5028e-03 1.6 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatAXPY               31 1.0 6.5468e-02 1.0 0.00e+00 0.0 3.6e+03 1.2e+02 5.0e+02  0  0  1  0  3   0  0  1  0  3     0
KSPGMRESOrthog      4889 1.0 7.6119e-01 3.9 2.15e+08 1.7 0.0e+00 0.0e+00 4.9e+03  2 17  0  0 33   2 17  0  0 33  2665
KSPSetUp              65 1.0 5.5337e-04 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPSolve              55 1.0 1.8826e+00 1.0 1.48e+09 2.5 2.9e+05 5.2e+02 1.0e+04  7 97 83 72 68   7 97 83 72 68  6018
PCSetUp               20 1.0 2.6713e-01 4.5 5.93e+07 5.3 0.0e+00 0.0e+00 5.4e+01  1  3  0  0  0   1  3  0  0  0  1275
PCSetUpOnBlocks       55 1.0 2.6609e-01 4.6 5.93e+07 5.3 0.0e+00 0.0e+00 5.0e+01  1  3  0  0  0   1  3  0  0  0  1280
PCApply             5128 1.0 1.0544e+00 2.3 1.01e+09 3.1 0.0e+00 0.0e+00 0.0e+00  3 60  0  0  0   3 60  0  0  0  6663
------------------------------------------------------------------------------------------------------------------------

Memory usage is given in bytes:

Object Type          Creations   Destructions     Memory  Descendants' Mem.
Reports information only for process 0.

--- Event Stage 0: Main Stage

              Vector   307            307      2686384     0
      Vector Scatter    47             47        48692     0
           Index Set   124            124       194436     0
   IS L to G Mapping     6              6         3384     0
              Matrix   133            133     21246352     0
       Krylov Solver     4              4        38720     0
      Preconditioner     4              4         3568     0
              Viewer     1              0            0     0
========================================================================================================================
Average time to get PetscTime(): 9.53674e-08
Average time for MPI_Barrier(): 5.19753e-06
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
| Time:           Thu Jan 31 22:20:06 2013                                                                             |
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
| libMesh Performance: Alive time=26.1485, Active time=25.7151                                                      |
 -------------------------------------------------------------------------------------------------------------------
| Event                                 nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                                 w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|-------------------------------------------------------------------------------------------------------------------|
|                                                                                                                   |
|                                                                                                                   |
| DofMap                                                                                                            |
|   add_neighbors_to_send_list()        2         0.1558      0.077902    0.3386      0.169294    0.61     1.32     |
|   build_constraint_matrix()           3816      0.0275      0.000007    0.0275      0.000007    0.11     0.11     |
|   build_sparsity()                    2         0.1346      0.067308    0.3274      0.163704    0.52     1.27     |
|   cnstrn_elem_mat_vec()               3816      0.0077      0.000002    0.0077      0.000002    0.03     0.03     |
|   create_dof_constraints()            2         0.1290      0.064492    0.6506      0.325277    0.50     2.53     |
|   distribute_dofs()                   2         0.1977      0.098869    0.7176      0.358781    0.77     2.79     |
|   dof_indices()                       69922     9.8695      0.000141    9.8695      0.000141    38.38    38.38    |
|   prepare_send_list()                 2         0.0028      0.001384    0.0028      0.001384    0.01     0.01     |
|   reinit()                            2         0.5018      0.250887    0.5018      0.250887    1.95     1.95     |
|                                                                                                                   |
| FE                                                                                                                |
|   compute_shape_functions()           36292     4.6595      0.000128    4.6595      0.000128    18.12    18.12    |
|   init_shape_functions()              3298      0.1944      0.000059    0.1944      0.000059    0.76     0.76     |
|   inverse_map()                       30558     0.7358      0.000024    0.7358      0.000024    2.86     2.86     |
|                                                                                                                   |
| FEMap                                                                                                             |
|   compute_affine_map()                36292     0.7330      0.000020    0.7330      0.000020    2.85     2.85     |
|   compute_face_map()                  3220      0.0620      0.000019    0.0620      0.000019    0.24     0.24     |
|   init_face_shape_functions()         20        0.0010      0.000051    0.0010      0.000051    0.00     0.00     |
|   init_reference_to_physical_map()    3298      0.1799      0.000055    0.1799      0.000055    0.70     0.70     |
|                                                                                                                   |
| Mesh                                                                                                              |
|   find_neighbors()                    1         0.1910      0.191009    0.1933      0.193296    0.74     0.75     |
|   renumber_nodes_and_elem()           2         0.0074      0.003703    0.0074      0.003703    0.03     0.03     |
|                                                                                                                   |
| MeshCommunication                                                                                                 |
|   assign_global_indices()             2         0.4642      0.232118    0.4723      0.236164    1.81     1.84     |
|   compute_hilbert_indices()           2         0.0882      0.044096    0.0882      0.044096    0.34     0.34     |
|   find_global_indices()               2         0.0334      0.016703    0.1317      0.065833    0.13     0.51     |
|   parallel_sort()                     2         0.0066      0.003286    0.0075      0.003743    0.03     0.03     |
|                                                                                                                   |
| MeshTools::Generation                                                                                             |
|   build_cube()                        1         0.0326      0.032621    0.0326      0.032621    0.13     0.13     |
|                                                                                                                   |
| MetisPartitioner                                                                                                  |
|   partition()                         1         0.3072      0.307205    0.3724      0.372434    1.19     1.45     |
|                                                                                                                   |
| Parallel                                                                                                          |
|   allgather()                         22        0.0123      0.000560    0.0124      0.000564    0.05     0.05     |
|   barrier()                           2         0.0001      0.000036    0.0001      0.000036    0.00     0.00     |
|   broadcast()                         45        0.0004      0.000010    0.0004      0.000010    0.00     0.00     |
|   max(bool)                           10        0.0001      0.000006    0.0001      0.000006    0.00     0.00     |
|   max(scalar)                         206       0.0069      0.000033    0.0069      0.000033    0.03     0.03     |
|   max(vector)                         45        0.0006      0.000014    0.0015      0.000034    0.00     0.01     |
|   maxloc(scalar)                      13        0.0517      0.003980    0.0517      0.003980    0.20     0.20     |
|   min(bool)                           220       0.0015      0.000007    0.0015      0.000007    0.01     0.01     |
|   min(scalar)                         180       0.0636      0.000354    0.0636      0.000354    0.25     0.25     |
|   min(vector)                         45        0.0008      0.000018    0.0024      0.000054    0.00     0.01     |
|   probe()                             478       0.0142      0.000030    0.0142      0.000030    0.06     0.06     |
|   receive()                           382       0.0038      0.000010    0.0161      0.000042    0.01     0.06     |
|   send()                              294       0.0011      0.000004    0.0011      0.000004    0.00     0.00     |
|   send_receive()                      298       0.0042      0.000014    0.0209      0.000070    0.02     0.08     |
|   sum()                               45        0.0023      0.000052    0.0028      0.000063    0.01     0.01     |
|                                                                                                                   |
| Parallel::Request                                                                                                 |
|   wait()                              294       0.0006      0.000002    0.0006      0.000002    0.00     0.00     |
|                                                                                                                   |
| Partitioner                                                                                                       |
|   set_node_processor_ids()            1         0.0111      0.011142    0.0199      0.019882    0.04     0.08     |
|   set_parent_processor_ids()          1         0.0087      0.008729    0.0087      0.008729    0.03     0.03     |
|                                                                                                                   |
| PetscLinearSolver                                                                                                 |
|   solve()                             55        2.0613      0.037479    2.0613      0.037479    8.02     8.02     |
|                                                                                                                   |
| PointLocatorTree                                                                                                  |
|   init(no master)                     1         0.0836      0.083619    0.0838      0.083842    0.33     0.33     |
|   operator()                          15        0.0036      0.000237    0.0040      0.000267    0.01     0.02     |
|                                                                                                                   |
| RBConstruction                                                                                                    |
|   add_scaled_matrix_and_vector()      10        1.5934      0.159336    6.0212      0.602120    6.20     23.42    |
|   clear()                             3         0.0024      0.000806    0.0024      0.000806    0.01     0.01     |
|   compute_Fq_representor_innerprods() 2         0.0077      0.003850    0.1038      0.051896    0.03     0.40     |
|   compute_max_error_bound()           9         0.0127      0.001414    0.5605      0.062279    0.05     2.18     |
|   enrich_RB_space()                   4         0.0017      0.000415    0.0017      0.000415    0.01     0.01     |
|   train_reduced_basis()               2         0.0048      0.002413    3.7099      1.854961    0.02     14.43    |
|   truth_assembly()                    4         0.0488      0.012196    0.0490      0.012256    0.19     0.19     |
|   truth_solve()                       4         0.0010      0.000253    0.2780      0.069495    0.00     1.08     |
|   update_RB_system_matrices()         8         0.0124      0.001554    0.0124      0.001554    0.05     0.05     |
|   update_residual_terms()             4         0.0595      0.014866    1.0399      0.259982    0.23     4.04     |
|                                                                                                                   |
| RBEIMConstruction                                                                                                 |
|   compute_best_fit_error()            100       0.0543      0.000543    0.0670      0.000670    0.21     0.26     |
|   enrich_RB_space()                   4         0.1100      0.027506    1.7005      0.425131    0.43     6.61     |
|   truth_solve()                       129       1.9898      0.015425    12.4257     0.096323    7.74     48.32    |
|   update_RB_system_matrices()         4         0.0012      0.000308    0.0092      0.002312    0.00     0.04     |
|                                                                                                                   |
| RBEIMEvaluation                                                                                                   |
|   rb_solve()                          426       0.0092      0.000022    0.0092      0.000022    0.04     0.04     |
|   write_offline_data_to_files()       1         0.0002      0.000202    0.0009      0.000889    0.00     0.00     |
|                                                                                                                   |
| RBEvaluation                                                                                                      |
|   clear()                             3         0.0001      0.000042    0.0001      0.000042    0.00     0.00     |
|   compute_residual_dual_norm()        420       0.4462      0.001062    0.4462      0.001062    1.74     1.74     |
|   rb_solve()                          420       0.0250      0.000059    0.4803      0.001144    0.10     1.87     |
|   resize_data_structures()            2         0.0008      0.000395    0.0008      0.000395    0.00     0.00     |
|   write_offline_data_to_files()       2         0.0014      0.000713    0.0014      0.000713    0.01     0.01     |
|   write_out_basis_functions()         2         0.0002      0.000089    0.7550      0.377493    0.00     2.94     |
|   write_out_vectors()                 2         0.2771      0.138529    0.7548      0.377403    1.08     2.94     |
 -------------------------------------------------------------------------------------------------------------------
| Totals:                               194774    25.7151                                         100.00            |
 -------------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example reduced_basis_ex6:
*  mpirun -np 12 example-devel -online_mode 0 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
***************************************************************
***************************************************************
* Running Example reduced_basis_ex6:
*  mpirun -np 12 example-devel -online_mode 1 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
***************************************************************
 
*** Warning, This code is untested, experimental, or likely to see future API changes: src/reduced_basis/rb_parametrized.C, line 41, compiled Jan 31 2013 at 21:51:32 ***
 EquationSystems
  n_systems()=2
   System #0, "EIM"
    Type "RBConstruction"
    Variables={ "x_comp_of_G" "y_comp_of_G" "z_comp_of_G" } 
    Finite Element Types="LAGRANGE", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" 
    Approximation Orders="FIRST", "THIRD" 
    n_dofs()=18513
    n_local_dofs()=1911
    n_constrained_dofs()=0
    n_local_constrained_dofs()=0
    n_vectors()=1
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 22.0431
      Average Off-Processor Bandwidth <= 3.38324
      Maximum  On-Processor Bandwidth <= 39
      Maximum Off-Processor Bandwidth <= 27
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
    n_local_dofs()=637
    n_constrained_dofs()=242
    n_local_constrained_dofs()=0
    n_vectors()=1
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 22.0431
      Average Off-Processor Bandwidth <= 3.38324
      Maximum  On-Processor Bandwidth <= 39
      Maximum Off-Processor Bandwidth <= 27
    DofMap Constraints
      Number of DoF Constraints = 242
      Average DoF Constraint Length= 0
      Number of Node Constraints = 0

 Mesh Information:
  mesh_dimension()=3
  spatial_dimension()=3
  n_nodes()=6171
    n_local_nodes()=637
  n_elem()=5000
    n_local_elem()=424
    n_active_elem()=5000
  n_subdomains()=1
  n_partitions()=12
  n_processors()=12
  n_threads()=1
  processor_id()=0

Bi: 5.000000e-03
curvature: 1.047200e+00
kappa: 1.300000e+00

************************************************************************************************************************
***             WIDEN YOUR WINDOW TO 120 CHARACTERS.  Use 'enscript -r -fCourier9' to print this document            ***
************************************************************************************************************************

---------------------------------------------- PETSc Performance Summary: ----------------------------------------------

/workspace/libmesh/examples/reduced_basis/reduced_basis_ex6/.libs/lt-example-devel on a intel-12. named hbar.ices.utexas.edu with 12 processors, by benkirk Thu Jan 31 22:20:10 2013
Using Petsc Release Version 3.3.0, Patch 2, Fri Jul 13 15:42:00 CDT 2012 

                         Max       Max/Min        Avg      Total 
Time (sec):           3.979e+00      1.00144   3.978e+00
Objects:              3.200e+01      1.00000   3.200e+01
Flops:                1.703e+04      1.71466   1.337e+04  1.604e+05
Flops/sec:            4.280e+03      1.71223   3.361e+03  4.033e+04
MPI Messages:         4.200e+01      3.50000   2.900e+01  3.480e+02
MPI Message Lengths:  2.426e+04      1.92694   6.824e+02  2.375e+05
MPI Reductions:       6.900e+01      1.00000

Flop counting convention: 1 flop = 1 real number operation of type (multiply/divide/add/subtract)
                            e.g., VecAXPY() for real vectors of length N --> 2N flops
                            and VecAXPY() for complex vectors of length N --> 8N flops

Summary of Stages:   ----- Time ------  ----- Flops -----  --- Messages ---  -- Message Lengths --  -- Reductions --
                        Avg     %Total     Avg     %Total   counts   %Total     Avg         %Total   counts   %Total 
 0:      Main Stage: 3.9778e+00 100.0%  1.6045e+05 100.0%  3.480e+02 100.0%  6.824e+02      100.0%  6.800e+01  98.6% 

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

VecCopy                2 1.0 1.7166e-05 1.7 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecSet                17 1.0 3.3379e-05 1.8 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecAXPY                7 1.0 9.0122e-05 1.6 1.70e+04 1.7 0.0e+00 0.0e+00 0.0e+00  0100  0  0  0   0100  0  0  0  1780
VecAssemblyBegin       9 1.0 4.8518e-04 1.5 0.00e+00 0.0 0.0e+00 0.0e+00 2.7e+01  0  0  0  0 39   0  0  0  0 40     0
VecAssemblyEnd         9 1.0 4.5061e-05 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecScatterBegin        2 1.0 1.0395e-04 1.2 0.00e+00 0.0 1.2e+02 1.4e+03 0.0e+00  0  0 33 67  0   0  0 33 67  0     0
VecScatterEnd          2 1.0 5.0068e-05 2.8 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatZeroEntries         4 1.0 2.5392e-04 2.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
------------------------------------------------------------------------------------------------------------------------

Memory usage is given in bytes:

Object Type          Creations   Destructions     Memory  Descendants' Mem.
Reports information only for process 0.

--- Event Stage 0: Main Stage

              Vector    17             17       182528     0
      Vector Scatter     2              2         2072     0
           Index Set     4              4         7680     0
   IS L to G Mapping     2              2         1128     0
              Matrix     6              6       922336     0
              Viewer     1              0            0     0
========================================================================================================================
Average time to get PetscTime(): 1.19209e-07
Average time for MPI_Barrier(): 4.40121e-05
Average time for zero size MPI_Send(): 3.76701e-05
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
| Time:           Thu Jan 31 22:20:10 2013                                                                             |
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
| libMesh Performance: Alive time=4.04264, Active time=3.84978                                                 |
 --------------------------------------------------------------------------------------------------------------
| Event                            nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                            w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|--------------------------------------------------------------------------------------------------------------|
|                                                                                                              |
|                                                                                                              |
| DofMap                                                                                                       |
|   add_neighbors_to_send_list()   2         0.1546      0.077317    0.3384      0.169207    4.02     8.79     |
|   build_sparsity()               2         0.1320      0.066003    0.3489      0.174450    3.43     9.06     |
|   create_dof_constraints()       2         0.1215      0.060731    0.6456      0.322815    3.16     16.77    |
|   distribute_dofs()              2         0.1949      0.097445    0.7207      0.360335    5.06     18.72    |
|   dof_indices()                  9290      1.0553      0.000114    1.0553      0.000114    27.41    27.41    |
|   prepare_send_list()            2         0.0027      0.001329    0.0027      0.001329    0.07     0.07     |
|   reinit()                       2         0.4938      0.246883    0.4938      0.246883    12.83    12.83    |
|                                                                                                              |
| EquationSystems                                                                                              |
|   build_solution_vector()        1         0.0191      0.019131    0.1989      0.198866    0.50     5.17     |
|                                                                                                              |
| ExodusII_IO                                                                                                  |
|   write_nodal_data()             1         0.0530      0.053012    0.0530      0.053012    1.38     1.38     |
|                                                                                                              |
| Mesh                                                                                                         |
|   find_neighbors()               1         0.1804      0.180357    0.3209      0.320882    4.68     8.34     |
|   renumber_nodes_and_elem()      2         0.0073      0.003647    0.0073      0.003647    0.19     0.19     |
|                                                                                                              |
| MeshCommunication                                                                                            |
|   assign_global_indices()        2         0.4587      0.229369    0.4743      0.237144    11.92    12.32    |
|   compute_hilbert_indices()      2         0.0883      0.044134    0.0883      0.044134    2.29     2.29     |
|   find_global_indices()          2         0.0338      0.016910    0.1314      0.065708    0.88     3.41     |
|   parallel_sort()                2         0.0066      0.003291    0.0078      0.003883    0.17     0.20     |
|                                                                                                              |
| MeshOutput                                                                                                   |
|   write_equation_systems()       1         0.0001      0.000147    0.2525      0.252530    0.00     6.56     |
|                                                                                                              |
| MeshTools::Generation                                                                                        |
|   build_cube()                   1         0.0610      0.061035    0.0610      0.061035    1.59     1.59     |
|                                                                                                              |
| MetisPartitioner                                                                                             |
|   partition()                    1         0.3017      0.301653    0.3665      0.366471    7.84     9.52     |
|                                                                                                              |
| Parallel                                                                                                     |
|   allgather()                    22        0.0317      0.001440    0.0318      0.001444    0.82     0.83     |
|   barrier()                      2         0.0000      0.000020    0.0000      0.000020    0.00     0.00     |
|   broadcast()                    48        0.0005      0.000010    0.0004      0.000008    0.01     0.01     |
|   max(bool)                      2         0.0001      0.000025    0.0001      0.000025    0.00     0.00     |
|   max(scalar)                    177       0.0089      0.000050    0.0089      0.000050    0.23     0.23     |
|   max(vector)                    39        0.0009      0.000023    0.0027      0.000069    0.02     0.07     |
|   min(bool)                      201       0.0031      0.000015    0.0031      0.000015    0.08     0.08     |
|   min(scalar)                    166       0.2373      0.001429    0.2373      0.001429    6.16     6.16     |
|   min(vector)                    39        0.0010      0.000026    0.0037      0.000095    0.03     0.10     |
|   probe()                        382       0.0108      0.000028    0.0108      0.000028    0.28     0.28     |
|   receive()                      338       0.0031      0.000009    0.0138      0.000041    0.08     0.36     |
|   send()                         338       0.0014      0.000004    0.0014      0.000004    0.04     0.04     |
|   send_receive()                 298       0.0041      0.000014    0.0192      0.000064    0.11     0.50     |
|   sum()                          48        0.0044      0.000093    0.0057      0.000119    0.12     0.15     |
|                                                                                                              |
| Parallel::Request                                                                                            |
|   wait()                         338       0.0007      0.000002    0.0007      0.000002    0.02     0.02     |
|                                                                                                              |
| Partitioner                                                                                                  |
|   set_node_processor_ids()       1         0.0111      0.011059    0.0273      0.027292    0.29     0.71     |
|   set_parent_processor_ids()     1         0.0092      0.009166    0.0092      0.009166    0.24     0.24     |
|                                                                                                              |
| RBConstruction                                                                                               |
|   clear()                        3         0.0011      0.000352    0.0011      0.000352    0.03     0.03     |
|   load_rb_solution()             2         0.0006      0.000295    0.0006      0.000295    0.02     0.02     |
|                                                                                                              |
| RBEIMEvaluation                                                                                              |
|   rb_solve()                     1         0.0135      0.013538    0.0135      0.013538    0.35     0.35     |
|   read_offline_data_from_files() 1         0.0002      0.000158    0.0011      0.001119    0.00     0.03     |
|                                                                                                              |
| RBEvaluation                                                                                                 |
|   clear()                        3         0.0002      0.000058    0.0002      0.000058    0.00     0.00     |
|   compute_residual_dual_norm()   1         0.0027      0.002679    0.0027      0.002679    0.07     0.07     |
|   rb_solve()                     1         0.0002      0.000166    0.0164      0.016384    0.00     0.43     |
|   read_in_basis_functions()      2         0.0002      0.000078    0.6132      0.306584    0.00     15.93    |
|   read_in_vectors()              2         0.1365      0.068235    0.6130      0.306504    3.54     15.92    |
|   read_offline_data_from_files() 2         0.0013      0.000645    0.0019      0.000955    0.03     0.05     |
|   resize_data_structures()       2         0.0006      0.000311    0.0006      0.000311    0.02     0.02     |
 --------------------------------------------------------------------------------------------------------------
| Totals:                          11778     3.8498                                          100.00            |
 --------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example reduced_basis_ex6:
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
