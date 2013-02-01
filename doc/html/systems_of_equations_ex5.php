<?php $root=""; ?>
<?php require($root."navigation.php"); ?>
<html>
<head>
  <?php load_style($root); ?>
</head>
 
<body>
 
<?php make_navigation("systems_of_equations_ex5",$root)?>
 
<div class="content">
<a name="comments"></a> 
<br><br><br> <h1> The source file systems_of_equations_ex5.C with comments: </h1> 
<div class = "comment">
<h1> Systems Example 5 - Linear Elastic Cantilever with Constraint </h1>
By David Knezevic

<br><br>In this example we extend systems_of_equations_ex4 to enforce a constraint.
We apply a uniform load on the top surface of the cantilever, and we
determine the traction on the right boundary in order to obtain zero
average vertical displacement on the right boundary of the domain.

<br><br>This constraint is enforced via a Lagrange multiplier (SCALAR variable).
The system we solve, therefore, is of the form:
a(u,v) + \lambda g(v) = f(v)
g(u) = 0
Here \lambda tells us the traction required to satisfy the constraint.


<br><br>C++ include files that we need
</div>

<div class ="fragment">
<pre>
        #include &lt;iostream&gt;
        #include &lt;algorithm&gt;
        #include &lt;math.h&gt;
        
</pre>
</div>
<div class = "comment">
libMesh includes
</div>

<div class ="fragment">
<pre>
        #include "libmesh/libmesh.h"
        #include "libmesh/mesh.h"
        #include "libmesh/mesh_generation.h"
        #include "libmesh/exodusII_io.h"
        #include "libmesh/linear_implicit_system.h"
        #include "libmesh/equation_systems.h"
        #include "libmesh/fe.h"
        #include "libmesh/quadrature_gauss.h"
        #include "libmesh/dof_map.h"
        #include "libmesh/sparse_matrix.h"
        #include "libmesh/numeric_vector.h"
        #include "libmesh/dense_matrix.h"
        #include "libmesh/dense_submatrix.h"
        #include "libmesh/dense_vector.h"
        #include "libmesh/dense_subvector.h"
        #include "libmesh/perf_log.h"
        #include "libmesh/elem.h"
        #include "libmesh/boundary_info.h"
        #include "libmesh/zero_function.h"
        #include "libmesh/dirichlet_boundaries.h"
        #include "libmesh/string_to_enum.h"
        #include "libmesh/getpot.h"
        
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
Matrix and right-hand side assemble
</div>

<div class ="fragment">
<pre>
        void assemble_elasticity(EquationSystems& es,
                                 const std::string& system_name);
        
</pre>
</div>
<div class = "comment">
Define the elasticity tensor, which is a fourth-order tensor
i.e. it has four indices i,j,k,l
</div>

<div class ="fragment">
<pre>
        Real eval_elasticity_tensor(unsigned int i,
                                    unsigned int j,
                                    unsigned int k,
                                    unsigned int l);
        
</pre>
</div>
<div class = "comment">
Begin the main program.
</div>

<div class ="fragment">
<pre>
        int main (int argc, char** argv)
        {
</pre>
</div>
<div class = "comment">
Initialize libMesh and any dependent libaries
</div>

<div class ="fragment">
<pre>
          LibMeshInit init (argc, argv);
        
</pre>
</div>
<div class = "comment">
Initialize the cantilever mesh
</div>

<div class ="fragment">
<pre>
          const unsigned int dim = 2;
        
</pre>
</div>
<div class = "comment">
Skip this 2D example if libMesh was compiled as 1D-only.
</div>

<div class ="fragment">
<pre>
          libmesh_example_assert(dim &lt;= LIBMESH_DIM, "2D support");
        
          Mesh mesh(dim);
          MeshTools::Generation::build_square (mesh,
                                               50, 10,
                                               0., 1.,
                                               0., 0.2,
                                               QUAD9);
        
          
</pre>
</div>
<div class = "comment">
Print information about the mesh to the screen.
</div>

<div class ="fragment">
<pre>
          mesh.print_info();
          
          
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
Declare the system and its variables.
Create a system named "Elasticity"
</div>

<div class ="fragment">
<pre>
          LinearImplicitSystem& system =
            equation_systems.add_system&lt;LinearImplicitSystem&gt; ("Elasticity");
        
          
</pre>
</div>
<div class = "comment">
Add two displacement variables, u and v, to the system
</div>

<div class ="fragment">
<pre>
          unsigned int u_var      = system.add_variable("u", SECOND, LAGRANGE);
          unsigned int v_var      = system.add_variable("v", SECOND, LAGRANGE);
          
</pre>
</div>
<div class = "comment">
Add a SCALAR variable for the Lagrange multiplier to enforce our constraint
</div>

<div class ="fragment">
<pre>
          unsigned int lambda_var = system.add_variable("lambda", FIRST, SCALAR);
        
        
          system.attach_assemble_function (assemble_elasticity);
        
</pre>
</div>
<div class = "comment">
Construct a Dirichlet boundary condition object
We impose a "clamped" boundary condition on the
"left" boundary, i.e. bc_id = 3
</div>

<div class ="fragment">
<pre>
          std::set&lt;boundary_id_type&gt; boundary_ids;
          boundary_ids.insert(3);
        
</pre>
</div>
<div class = "comment">
Create a vector storing the variable numbers which the BC applies to
</div>

<div class ="fragment">
<pre>
          std::vector&lt;unsigned int&gt; variables(2);
          variables[0] = u_var; variables[1] = v_var;
          
</pre>
</div>
<div class = "comment">
Create a ZeroFunction to initialize dirichlet_bc
</div>

<div class ="fragment">
<pre>
          ZeroFunction&lt;&gt; zf;
          
          DirichletBoundary dirichlet_bc(boundary_ids,
                                         variables,
                                         &zf);
        
</pre>
</div>
<div class = "comment">
We must add the Dirichlet boundary condition _before_ 
we call equation_systems.init()
</div>

<div class ="fragment">
<pre>
          system.get_dof_map().add_dirichlet_boundary(dirichlet_bc);
          
</pre>
</div>
<div class = "comment">
Initialize the data structures for the equation system.
</div>

<div class ="fragment">
<pre>
          equation_systems.init();
        
</pre>
</div>
<div class = "comment">
Print information about the system to the screen.
</div>

<div class ="fragment">
<pre>
          equation_systems.print_info();
        
</pre>
</div>
<div class = "comment">
Solve the system
</div>

<div class ="fragment">
<pre>
          system.solve();
        
</pre>
</div>
<div class = "comment">
Plot the solution
</div>

<div class ="fragment">
<pre>
        #ifdef LIBMESH_HAVE_EXODUS_API
          ExodusII_IO (mesh).write_equation_systems("displacement.e",equation_systems);
        #endif // #ifdef LIBMESH_HAVE_EXODUS_API
        
</pre>
</div>
<div class = "comment">
All done.  
</div>

<div class ="fragment">
<pre>
          return 0;
        }
        
        
        void assemble_elasticity(EquationSystems& es,
                                 const std::string& system_name)
        {
          libmesh_assert_equal_to (system_name, "Elasticity");
          
          const MeshBase& mesh = es.get_mesh();
        
          const unsigned int dim = mesh.mesh_dimension();
        
          LinearImplicitSystem& system = es.get_system&lt;LinearImplicitSystem&gt;("Elasticity");
        
          const unsigned int u_var = system.variable_number ("u");
          const unsigned int v_var = system.variable_number ("v");
          const unsigned int lambda_var = system.variable_number ("lambda");
        
          const DofMap& dof_map = system.get_dof_map();
          FEType fe_type = dof_map.variable_type(0);
          AutoPtr&lt;FEBase&gt; fe (FEBase::build(dim, fe_type));
          QGauss qrule (dim, fe_type.default_quadrature_order());
          fe-&gt;attach_quadrature_rule (&qrule);
        
          AutoPtr&lt;FEBase&gt; fe_face (FEBase::build(dim, fe_type));
          QGauss qface(dim-1, fe_type.default_quadrature_order());
          fe_face-&gt;attach_quadrature_rule (&qface);
        
          const std::vector&lt;Real&gt;& JxW = fe-&gt;get_JxW();
          const std::vector&lt;std::vector&lt;RealGradient&gt; &gt;& dphi = fe-&gt;get_dphi();
        
          DenseMatrix&lt;Number&gt; Ke;
          DenseVector&lt;Number&gt; Fe;
        
          DenseSubMatrix&lt;Number&gt;
            Kuu(Ke), Kuv(Ke),
            Kvu(Ke), Kvv(Ke);
          DenseSubMatrix&lt;Number&gt; Klambda_v(Ke), Kv_lambda(Ke);
        
          DenseSubVector&lt;Number&gt;
            Fu(Fe),
            Fv(Fe);
        
          std::vector&lt;dof_id_type&gt; dof_indices;
          std::vector&lt;dof_id_type&gt; dof_indices_u;
          std::vector&lt;dof_id_type&gt; dof_indices_v;
          std::vector&lt;dof_id_type&gt; dof_indices_lambda;
        
          MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
          const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
        
          for ( ; el != end_el; ++el)
            {
              const Elem* elem = *el;
        
              dof_map.dof_indices (elem, dof_indices);
              dof_map.dof_indices (elem, dof_indices_u, u_var);
              dof_map.dof_indices (elem, dof_indices_v, v_var);
              dof_map.dof_indices (elem, dof_indices_lambda, lambda_var);
        
              const unsigned int n_dofs   = dof_indices.size();
              const unsigned int n_u_dofs = dof_indices_u.size(); 
              const unsigned int n_v_dofs = dof_indices_v.size();
              const unsigned int n_lambda_dofs = dof_indices_lambda.size();
        
              fe-&gt;reinit (elem);
        
              Ke.resize (n_dofs, n_dofs);
              Fe.resize (n_dofs);
        
              Kuu.reposition (u_var*n_u_dofs, u_var*n_u_dofs, n_u_dofs, n_u_dofs);
              Kuv.reposition (u_var*n_u_dofs, v_var*n_u_dofs, n_u_dofs, n_v_dofs);
              
              Kvu.reposition (v_var*n_v_dofs, u_var*n_v_dofs, n_v_dofs, n_u_dofs);
              Kvv.reposition (v_var*n_v_dofs, v_var*n_v_dofs, n_v_dofs, n_v_dofs);
        
</pre>
</div>
<div class = "comment">
Also, add a row and a column to enforce the constraint
</div>

<div class ="fragment">
<pre>
              Kv_lambda.reposition (v_var*n_u_dofs, v_var*n_u_dofs+n_v_dofs, n_v_dofs, 1);
              Klambda_v.reposition (v_var*n_v_dofs+n_v_dofs, v_var*n_v_dofs, 1, n_v_dofs);
        
              Fu.reposition (u_var*n_u_dofs, n_u_dofs);
              Fv.reposition (v_var*n_u_dofs, n_v_dofs);
        
              for (unsigned int qp=0; qp&lt;qrule.n_points(); qp++)
              {
                  for (unsigned int i=0; i&lt;n_u_dofs; i++)
                    for (unsigned int j=0; j&lt;n_u_dofs; j++)
                    {
</pre>
</div>
<div class = "comment">
Tensor indices
</div>

<div class ="fragment">
<pre>
                      unsigned int C_i, C_j, C_k, C_l;
                      C_i=0, C_k=0;
        
        
                      C_j=0, C_l=0;
                      Kuu(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i,C_j,C_k,C_l) * dphi[i][qp](C_j)*dphi[j][qp](C_l));
                      
                      C_j=1, C_l=0;
                      Kuu(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i,C_j,C_k,C_l) * dphi[i][qp](C_j)*dphi[j][qp](C_l));
        
                      C_j=0, C_l=1;
                      Kuu(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i,C_j,C_k,C_l) * dphi[i][qp](C_j)*dphi[j][qp](C_l));
        
                      C_j=1, C_l=1;
                      Kuu(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i,C_j,C_k,C_l) * dphi[i][qp](C_j)*dphi[j][qp](C_l));
                    }
        
                  for (unsigned int i=0; i&lt;n_u_dofs; i++)
                    for (unsigned int j=0; j&lt;n_v_dofs; j++)
                    {
</pre>
</div>
<div class = "comment">
Tensor indices
</div>

<div class ="fragment">
<pre>
                      unsigned int C_i, C_j, C_k, C_l;
                      C_i=0, C_k=1;
        
        
                      C_j=0, C_l=0;
                      Kuv(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i,C_j,C_k,C_l) * dphi[i][qp](C_j)*dphi[j][qp](C_l));
                      
                      C_j=1, C_l=0;
                      Kuv(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i,C_j,C_k,C_l) * dphi[i][qp](C_j)*dphi[j][qp](C_l));
        
                      C_j=0, C_l=1;
                      Kuv(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i,C_j,C_k,C_l) * dphi[i][qp](C_j)*dphi[j][qp](C_l));
        
                      C_j=1, C_l=1;
                      Kuv(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i,C_j,C_k,C_l) * dphi[i][qp](C_j)*dphi[j][qp](C_l));
                    }
        
                  for (unsigned int i=0; i&lt;n_v_dofs; i++)
                    for (unsigned int j=0; j&lt;n_u_dofs; j++)
                    {
</pre>
</div>
<div class = "comment">
Tensor indices
</div>

<div class ="fragment">
<pre>
                      unsigned int C_i, C_j, C_k, C_l;
                      C_i=1, C_k=0;
        
        
                      C_j=0, C_l=0;
                      Kvu(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i,C_j,C_k,C_l) * dphi[i][qp](C_j)*dphi[j][qp](C_l));
                      
                      C_j=1, C_l=0;
                      Kvu(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i,C_j,C_k,C_l) * dphi[i][qp](C_j)*dphi[j][qp](C_l));
        
                      C_j=0, C_l=1;
                      Kvu(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i,C_j,C_k,C_l) * dphi[i][qp](C_j)*dphi[j][qp](C_l));
        
                      C_j=1, C_l=1;
                      Kvu(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i,C_j,C_k,C_l) * dphi[i][qp](C_j)*dphi[j][qp](C_l));
                    }
        
                  for (unsigned int i=0; i&lt;n_v_dofs; i++)
                    for (unsigned int j=0; j&lt;n_v_dofs; j++)
                    {
</pre>
</div>
<div class = "comment">
Tensor indices
</div>

<div class ="fragment">
<pre>
                      unsigned int C_i, C_j, C_k, C_l;
                      C_i=1, C_k=1;
        
        
                      C_j=0, C_l=0;
                      Kvv(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i,C_j,C_k,C_l) * dphi[i][qp](C_j)*dphi[j][qp](C_l));
                      
                      C_j=1, C_l=0;
                      Kvv(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i,C_j,C_k,C_l) * dphi[i][qp](C_j)*dphi[j][qp](C_l));
        
                      C_j=0, C_l=1;
                      Kvv(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i,C_j,C_k,C_l) * dphi[i][qp](C_j)*dphi[j][qp](C_l));
        
                      C_j=1, C_l=1;
                      Kvv(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i,C_j,C_k,C_l) * dphi[i][qp](C_j)*dphi[j][qp](C_l));
                    }
              }
        
              {
                for (unsigned int side=0; side&lt;elem-&gt;n_sides(); side++)
                  if (elem-&gt;neighbor(side) == NULL)
                    {
                      const std::vector&lt;boundary_id_type&gt; bc_ids =
                        mesh.boundary_info-&gt;boundary_ids (elem,side);
        
                      const std::vector&lt;std::vector&lt;Real&gt; &gt;&  phi_face = fe_face-&gt;get_phi();
                      const std::vector&lt;Real&gt;& JxW_face = fe_face-&gt;get_JxW();
        
                      fe_face-&gt;reinit(elem, side);
        
                      for (std::vector&lt;boundary_id_type&gt;::const_iterator b =
                           bc_ids.begin(); b != bc_ids.end(); ++b)
                        {
                          const boundary_id_type bc_id = *b;
                          for (unsigned int qp=0; qp&lt;qface.n_points(); qp++)
                            {
</pre>
</div>
<div class = "comment">
Add the loading
</div>

<div class ="fragment">
<pre>
                              if( bc_id == 2 )
                                {
                                  for (unsigned int i=0; i&lt;n_v_dofs; i++)
                                    {
                                      Fv(i) += JxW_face[qp]* (-1.) * phi_face[i][qp];
                                    }
                                }
        
</pre>
</div>
<div class = "comment">
Add the constraint contributions
</div>

<div class ="fragment">
<pre>
                              if( bc_id == 1 )
                                {
                                  for (unsigned int i=0; i&lt;n_v_dofs; i++)
                                    for (unsigned int j=0; j&lt;n_lambda_dofs; j++)
                                    {
                                      Kv_lambda(i,j) += JxW_face[qp]* (-1.) * phi_face[i][qp];
                                    }
                                
                                  for (unsigned int i=0; i&lt;n_lambda_dofs; i++)
                                    for (unsigned int j=0; j&lt;n_v_dofs; j++)
                                      {
                                        Klambda_v(i,j) += JxW_face[qp]* (-1.) * phi_face[j][qp];
                                      }
                                }
                            }
                        }
                    }
              } 
        
              dof_map.constrain_element_matrix_and_vector (Ke, Fe, dof_indices);
              
              system.matrix-&gt;add_matrix (Ke, dof_indices);
              system.rhs-&gt;add_vector    (Fe, dof_indices);
            }
        }
        
        Real eval_elasticity_tensor(unsigned int i,
                                    unsigned int j,
                                    unsigned int k,
                                    unsigned int l)
        {
</pre>
</div>
<div class = "comment">
Define the Poisson ratio
</div>

<div class ="fragment">
<pre>
          const Real nu = 0.3;
          
</pre>
</div>
<div class = "comment">
Define the Lame constants (lambda_1 and lambda_2) based on Poisson ratio
</div>

<div class ="fragment">
<pre>
          const Real lambda_1 = nu / ( (1. + nu) * (1. - 2.*nu) );
          const Real lambda_2 = 0.5 / (1 + nu);
        
</pre>
</div>
<div class = "comment">
Define the Kronecker delta functions that we need here
</div>

<div class ="fragment">
<pre>
          Real delta_ij = (i == j) ? 1. : 0.;
          Real delta_il = (i == l) ? 1. : 0.;
          Real delta_ik = (i == k) ? 1. : 0.;
          Real delta_jl = (j == l) ? 1. : 0.;
          Real delta_jk = (j == k) ? 1. : 0.;
          Real delta_kl = (k == l) ? 1. : 0.;
          
          return lambda_1 * delta_ij * delta_kl + lambda_2 * (delta_ik * delta_jl + delta_il * delta_jk);
        }
</pre>
</div>

<a name="nocomments"></a> 
<br><br><br> <h1> The source file systems_of_equations_ex5.C without comments: </h1> 
<pre> 
  
  #include &lt;iostream&gt;
  #include &lt;algorithm&gt;
  #include &lt;math.h&gt;
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/libmesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/mesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/mesh_generation.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/exodusII_io.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/linear_implicit_system.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/equation_systems.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/fe.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/quadrature_gauss.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/dof_map.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/sparse_matrix.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/numeric_vector.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/dense_matrix.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/dense_submatrix.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/dense_vector.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/dense_subvector.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/perf_log.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/elem.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/boundary_info.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/zero_function.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/dirichlet_boundaries.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/string_to_enum.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/getpot.h&quot;</FONT></B>
  
  using namespace libMesh;
  
  <B><FONT COLOR="#228B22">void</FONT></B> assemble_elasticity(EquationSystems&amp; es,
                           <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp; system_name);
  
  Real eval_elasticity_tensor(<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i,
                              <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j,
                              <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> k,
                              <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> l);
  
  <B><FONT COLOR="#228B22">int</FONT></B> main (<B><FONT COLOR="#228B22">int</FONT></B> argc, <B><FONT COLOR="#228B22">char</FONT></B>** argv)
  {
    LibMeshInit init (argc, argv);
  
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> dim = 2;
  
    libmesh_example_assert(dim &lt;= LIBMESH_DIM, <B><FONT COLOR="#BC8F8F">&quot;2D support&quot;</FONT></B>);
  
    Mesh mesh(dim);
    <B><FONT COLOR="#5F9EA0">MeshTools</FONT></B>::Generation::build_square (mesh,
                                         50, 10,
                                         0., 1.,
                                         0., 0.2,
                                         QUAD9);
  
    
    mesh.print_info();
    
    
    EquationSystems equation_systems (mesh);
    
    LinearImplicitSystem&amp; system =
      equation_systems.add_system&lt;LinearImplicitSystem&gt; (<B><FONT COLOR="#BC8F8F">&quot;Elasticity&quot;</FONT></B>);
  
    
    <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> u_var      = system.add_variable(<B><FONT COLOR="#BC8F8F">&quot;u&quot;</FONT></B>, SECOND, LAGRANGE);
    <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> v_var      = system.add_variable(<B><FONT COLOR="#BC8F8F">&quot;v&quot;</FONT></B>, SECOND, LAGRANGE);
    
    <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> lambda_var = system.add_variable(<B><FONT COLOR="#BC8F8F">&quot;lambda&quot;</FONT></B>, FIRST, SCALAR);
  
  
    system.attach_assemble_function (assemble_elasticity);
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::set&lt;boundary_id_type&gt; boundary_ids;
    boundary_ids.insert(3);
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B>&gt; variables(2);
    variables[0] = u_var; variables[1] = v_var;
    
    ZeroFunction&lt;&gt; zf;
    
    DirichletBoundary dirichlet_bc(boundary_ids,
                                   variables,
                                   &amp;zf);
  
    system.get_dof_map().add_dirichlet_boundary(dirichlet_bc);
    
    equation_systems.init();
  
    equation_systems.print_info();
  
    system.solve();
  
  #ifdef LIBMESH_HAVE_EXODUS_API
    ExodusII_IO (mesh).write_equation_systems(<B><FONT COLOR="#BC8F8F">&quot;displacement.e&quot;</FONT></B>,equation_systems);
  #endif <I><FONT COLOR="#B22222">// #ifdef LIBMESH_HAVE_EXODUS_API
</FONT></I>  
    <B><FONT COLOR="#A020F0">return</FONT></B> 0;
  }
  
  
  <B><FONT COLOR="#228B22">void</FONT></B> assemble_elasticity(EquationSystems&amp; es,
                           <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp; system_name)
  {
    libmesh_assert_equal_to (system_name, <B><FONT COLOR="#BC8F8F">&quot;Elasticity&quot;</FONT></B>);
    
    <B><FONT COLOR="#228B22">const</FONT></B> MeshBase&amp; mesh = es.get_mesh();
  
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> dim = mesh.mesh_dimension();
  
    LinearImplicitSystem&amp; system = es.get_system&lt;LinearImplicitSystem&gt;(<B><FONT COLOR="#BC8F8F">&quot;Elasticity&quot;</FONT></B>);
  
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> u_var = system.variable_number (<B><FONT COLOR="#BC8F8F">&quot;u&quot;</FONT></B>);
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> v_var = system.variable_number (<B><FONT COLOR="#BC8F8F">&quot;v&quot;</FONT></B>);
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> lambda_var = system.variable_number (<B><FONT COLOR="#BC8F8F">&quot;lambda&quot;</FONT></B>);
  
    <B><FONT COLOR="#228B22">const</FONT></B> DofMap&amp; dof_map = system.get_dof_map();
    FEType fe_type = dof_map.variable_type(0);
    AutoPtr&lt;FEBase&gt; fe (FEBase::build(dim, fe_type));
    QGauss qrule (dim, fe_type.default_quadrature_order());
    fe-&gt;attach_quadrature_rule (&amp;qrule);
  
    AutoPtr&lt;FEBase&gt; fe_face (FEBase::build(dim, fe_type));
    QGauss qface(dim-1, fe_type.default_quadrature_order());
    fe_face-&gt;attach_quadrature_rule (&amp;qface);
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Real&gt;&amp; JxW = fe-&gt;get_JxW();
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;RealGradient&gt; &gt;&amp; dphi = fe-&gt;get_dphi();
  
    DenseMatrix&lt;Number&gt; Ke;
    DenseVector&lt;Number&gt; Fe;
  
    DenseSubMatrix&lt;Number&gt;
      Kuu(Ke), Kuv(Ke),
      Kvu(Ke), Kvv(Ke);
    DenseSubMatrix&lt;Number&gt; Klambda_v(Ke), Kv_lambda(Ke);
  
    DenseSubVector&lt;Number&gt;
      Fu(Fe),
      Fv(Fe);
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;dof_id_type&gt; dof_indices;
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;dof_id_type&gt; dof_indices_u;
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;dof_id_type&gt; dof_indices_v;
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;dof_id_type&gt; dof_indices_lambda;
  
    <B><FONT COLOR="#5F9EA0">MeshBase</FONT></B>::const_element_iterator       el     = mesh.active_local_elements_begin();
    <B><FONT COLOR="#228B22">const</FONT></B> MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
  
    <B><FONT COLOR="#A020F0">for</FONT></B> ( ; el != end_el; ++el)
      {
        <B><FONT COLOR="#228B22">const</FONT></B> Elem* elem = *el;
  
        dof_map.dof_indices (elem, dof_indices);
        dof_map.dof_indices (elem, dof_indices_u, u_var);
        dof_map.dof_indices (elem, dof_indices_v, v_var);
        dof_map.dof_indices (elem, dof_indices_lambda, lambda_var);
  
        <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_dofs   = dof_indices.size();
        <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_u_dofs = dof_indices_u.size(); 
        <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_v_dofs = dof_indices_v.size();
        <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_lambda_dofs = dof_indices_lambda.size();
  
        fe-&gt;reinit (elem);
  
        Ke.resize (n_dofs, n_dofs);
        Fe.resize (n_dofs);
  
        Kuu.reposition (u_var*n_u_dofs, u_var*n_u_dofs, n_u_dofs, n_u_dofs);
        Kuv.reposition (u_var*n_u_dofs, v_var*n_u_dofs, n_u_dofs, n_v_dofs);
        
        Kvu.reposition (v_var*n_v_dofs, u_var*n_v_dofs, n_v_dofs, n_u_dofs);
        Kvv.reposition (v_var*n_v_dofs, v_var*n_v_dofs, n_v_dofs, n_v_dofs);
  
        Kv_lambda.reposition (v_var*n_u_dofs, v_var*n_u_dofs+n_v_dofs, n_v_dofs, 1);
        Klambda_v.reposition (v_var*n_v_dofs+n_v_dofs, v_var*n_v_dofs, 1, n_v_dofs);
  
        Fu.reposition (u_var*n_u_dofs, n_u_dofs);
        Fv.reposition (v_var*n_u_dofs, n_v_dofs);
  
        <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> qp=0; qp&lt;qrule.n_points(); qp++)
        {
            <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i&lt;n_u_dofs; i++)
              <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j=0; j&lt;n_u_dofs; j++)
              {
                <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> C_i, C_j, C_k, C_l;
                C_i=0, C_k=0;
  
  
                C_j=0, C_l=0;
                Kuu(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i,C_j,C_k,C_l) * dphi[i][qp](C_j)*dphi[j][qp](C_l));
                
                C_j=1, C_l=0;
                Kuu(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i,C_j,C_k,C_l) * dphi[i][qp](C_j)*dphi[j][qp](C_l));
  
                C_j=0, C_l=1;
                Kuu(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i,C_j,C_k,C_l) * dphi[i][qp](C_j)*dphi[j][qp](C_l));
  
                C_j=1, C_l=1;
                Kuu(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i,C_j,C_k,C_l) * dphi[i][qp](C_j)*dphi[j][qp](C_l));
              }
  
            <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i&lt;n_u_dofs; i++)
              <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j=0; j&lt;n_v_dofs; j++)
              {
                <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> C_i, C_j, C_k, C_l;
                C_i=0, C_k=1;
  
  
                C_j=0, C_l=0;
                Kuv(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i,C_j,C_k,C_l) * dphi[i][qp](C_j)*dphi[j][qp](C_l));
                
                C_j=1, C_l=0;
                Kuv(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i,C_j,C_k,C_l) * dphi[i][qp](C_j)*dphi[j][qp](C_l));
  
                C_j=0, C_l=1;
                Kuv(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i,C_j,C_k,C_l) * dphi[i][qp](C_j)*dphi[j][qp](C_l));
  
                C_j=1, C_l=1;
                Kuv(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i,C_j,C_k,C_l) * dphi[i][qp](C_j)*dphi[j][qp](C_l));
              }
  
            <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i&lt;n_v_dofs; i++)
              <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j=0; j&lt;n_u_dofs; j++)
              {
                <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> C_i, C_j, C_k, C_l;
                C_i=1, C_k=0;
  
  
                C_j=0, C_l=0;
                Kvu(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i,C_j,C_k,C_l) * dphi[i][qp](C_j)*dphi[j][qp](C_l));
                
                C_j=1, C_l=0;
                Kvu(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i,C_j,C_k,C_l) * dphi[i][qp](C_j)*dphi[j][qp](C_l));
  
                C_j=0, C_l=1;
                Kvu(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i,C_j,C_k,C_l) * dphi[i][qp](C_j)*dphi[j][qp](C_l));
  
                C_j=1, C_l=1;
                Kvu(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i,C_j,C_k,C_l) * dphi[i][qp](C_j)*dphi[j][qp](C_l));
              }
  
            <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i&lt;n_v_dofs; i++)
              <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j=0; j&lt;n_v_dofs; j++)
              {
                <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> C_i, C_j, C_k, C_l;
                C_i=1, C_k=1;
  
  
                C_j=0, C_l=0;
                Kvv(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i,C_j,C_k,C_l) * dphi[i][qp](C_j)*dphi[j][qp](C_l));
                
                C_j=1, C_l=0;
                Kvv(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i,C_j,C_k,C_l) * dphi[i][qp](C_j)*dphi[j][qp](C_l));
  
                C_j=0, C_l=1;
                Kvv(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i,C_j,C_k,C_l) * dphi[i][qp](C_j)*dphi[j][qp](C_l));
  
                C_j=1, C_l=1;
                Kvv(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i,C_j,C_k,C_l) * dphi[i][qp](C_j)*dphi[j][qp](C_l));
              }
        }
  
        {
          <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> side=0; side&lt;elem-&gt;n_sides(); side++)
            <B><FONT COLOR="#A020F0">if</FONT></B> (elem-&gt;neighbor(side) == NULL)
              {
                <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;boundary_id_type&gt; bc_ids =
                  mesh.boundary_info-&gt;boundary_ids (elem,side);
  
                <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;Real&gt; &gt;&amp;  phi_face = fe_face-&gt;get_phi();
                <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Real&gt;&amp; JxW_face = fe_face-&gt;get_JxW();
  
                fe_face-&gt;reinit(elem, side);
  
                <B><FONT COLOR="#A020F0">for</FONT></B> (std::vector&lt;boundary_id_type&gt;::const_iterator b =
                     bc_ids.begin(); b != bc_ids.end(); ++b)
                  {
                    <B><FONT COLOR="#228B22">const</FONT></B> boundary_id_type bc_id = *b;
                    <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> qp=0; qp&lt;qface.n_points(); qp++)
                      {
                        <B><FONT COLOR="#A020F0">if</FONT></B>( bc_id == 2 )
                          {
                            <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i&lt;n_v_dofs; i++)
                              {
                                Fv(i) += JxW_face[qp]* (-1.) * phi_face[i][qp];
                              }
                          }
  
                        <B><FONT COLOR="#A020F0">if</FONT></B>( bc_id == 1 )
                          {
                            <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i&lt;n_v_dofs; i++)
                              <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j=0; j&lt;n_lambda_dofs; j++)
                              {
                                Kv_lambda(i,j) += JxW_face[qp]* (-1.) * phi_face[i][qp];
                              }
                          
                            <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i&lt;n_lambda_dofs; i++)
                              <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j=0; j&lt;n_v_dofs; j++)
                                {
                                  Klambda_v(i,j) += JxW_face[qp]* (-1.) * phi_face[j][qp];
                                }
                          }
                      }
                  }
              }
        } 
  
        dof_map.constrain_element_matrix_and_vector (Ke, Fe, dof_indices);
        
        system.matrix-&gt;add_matrix (Ke, dof_indices);
        system.rhs-&gt;add_vector    (Fe, dof_indices);
      }
  }
  
  Real eval_elasticity_tensor(<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i,
                              <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j,
                              <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> k,
                              <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> l)
  {
    <B><FONT COLOR="#228B22">const</FONT></B> Real nu = 0.3;
    
    <B><FONT COLOR="#228B22">const</FONT></B> Real lambda_1 = nu / ( (1. + nu) * (1. - 2.*nu) );
    <B><FONT COLOR="#228B22">const</FONT></B> Real lambda_2 = 0.5 / (1 + nu);
  
    Real delta_ij = (i == j) ? 1. : 0.;
    Real delta_il = (i == l) ? 1. : 0.;
    Real delta_ik = (i == k) ? 1. : 0.;
    Real delta_jl = (j == l) ? 1. : 0.;
    Real delta_jk = (j == k) ? 1. : 0.;
    Real delta_kl = (k == l) ? 1. : 0.;
    
    <B><FONT COLOR="#A020F0">return</FONT></B> lambda_1 * delta_ij * delta_kl + lambda_2 * (delta_ik * delta_jl + delta_il * delta_jk);
  }
</pre> 
<a name="output"></a> 
<br><br><br> <h1> The console output of the program: </h1> 
<pre>
***************************************************************
* Running Example systems_of_equations_ex5:
*  mpirun -np 12 example-devel -ksp_type cg -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
***************************************************************
 
 Mesh Information:
  mesh_dimension()=2
  spatial_dimension()=3
  n_nodes()=2121
    n_local_nodes()=199
  n_elem()=500
    n_local_elem()=42
    n_active_elem()=500
  n_subdomains()=1
  n_partitions()=12
  n_processors()=12
  n_threads()=1
  processor_id()=0

 EquationSystems
  n_systems()=1
   System #0, "Elasticity"
    Type "LinearImplicit"
    Variables={ "u" "v" } "lambda" 
    Finite Element Types="LAGRANGE", "JACOBI_20_00" "SCALAR", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" "CARTESIAN" 
    Approximation Orders="SECOND", "THIRD" "FIRST", "THIRD" 
    n_dofs()=4243
    n_local_dofs()=398
    n_constrained_dofs()=42
    n_local_constrained_dofs()=0
    n_vectors()=1
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 28.3311
      Average Off-Processor Bandwidth <= 5.28918
      Maximum  On-Processor Bandwidth <= 329
      Maximum Off-Processor Bandwidth <= 3914
    DofMap Constraints
      Number of DoF Constraints = 42
      Average DoF Constraint Length= 0
      Number of Node Constraints = 0

************************************************************************************************************************
***             WIDEN YOUR WINDOW TO 120 CHARACTERS.  Use 'enscript -r -fCourier9' to print this document            ***
************************************************************************************************************************

---------------------------------------------- PETSc Performance Summary: ----------------------------------------------

/workspace/libmesh/examples/systems_of_equations/systems_of_equations_ex5/.libs/lt-example-devel on a intel-12. named hbar.ices.utexas.edu with 12 processors, by benkirk Thu Jan 31 22:15:13 2013
Using Petsc Release Version 3.3.0, Patch 2, Fri Jul 13 15:42:00 CDT 2012 

                         Max       Max/Min        Avg      Total 
Time (sec):           5.721e-01      1.01254   5.709e-01
Objects:              3.300e+01      1.00000   3.300e+01
Flops:                3.058e+06      2.95343   2.170e+06  2.604e+07
Flops/sec:            5.344e+06      2.95342   3.803e+06  4.563e+07
MPI Messages:         1.110e+02      5.69231   4.300e+01  5.160e+02
MPI Message Lengths:  1.487e+05      5.72216   1.058e+03  5.457e+05
MPI Reductions:       5.600e+01      1.00000

Flop counting convention: 1 flop = 1 real number operation of type (multiply/divide/add/subtract)
                            e.g., VecAXPY() for real vectors of length N --> 2N flops
                            and VecAXPY() for complex vectors of length N --> 8N flops

Summary of Stages:   ----- Time ------  ----- Flops -----  --- Messages ---  -- Message Lengths --  -- Reductions --
                        Avg     %Total     Avg     %Total   counts   %Total     Avg         %Total   counts   %Total 
 0:      Main Stage: 5.7089e-01 100.0%  2.6045e+07 100.0%  5.160e+02 100.0%  1.058e+03      100.0%  5.500e+01  98.2% 

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

VecTDot                5 1.0 1.0109e-04 2.1 3.98e+03 1.3 0.0e+00 0.0e+00 5.0e+00  0  0  0  0  9   0  0  0  0  9   419
VecNorm                4 1.0 1.0784e-02 1.2 3.18e+03 1.3 0.0e+00 0.0e+00 4.0e+00  2  0  0  0  7   2  0  0  0  7     3
VecCopy                2 1.0 4.0531e-06 1.3 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecSet                 9 1.0 1.0729e-05 2.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecAXPY                4 1.0 5.4121e-05 1.8 3.18e+03 1.3 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0   627
VecAYPX                2 1.0 5.2452e-06 1.4 1.19e+03 1.3 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0  2427
VecAssemblyBegin       3 1.0 1.9360e-04 1.1 0.00e+00 0.0 5.6e+01 2.7e+02 9.0e+00  0  0 11  3 16   0  0 11  3 16     0
VecAssemblyEnd         3 1.0 6.2227e-05 1.9 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecScatterBegin        4 1.0 8.1062e-05 1.4 0.00e+00 0.0 1.9e+02 7.8e+02 0.0e+00  0  0 37 27  0   0  0 37 27  0     0
VecScatterEnd          4 1.0 4.4339e-03315.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatMult                3 1.0 4.6039e-0325.9 8.38e+04 1.5 1.5e+02 8.6e+02 0.0e+00  0  3 29 24  0   0  3 29 24  0   178
MatSolve               4 1.0 3.0780e-04 1.7 3.36e+05 2.0 0.0e+00 0.0e+00 0.0e+00  0 13  0  0  0   0 13  0  0  0 10633
MatLUFactorNum         1 1.0 2.5170e-03 2.9 2.64e+06 3.3 0.0e+00 0.0e+00 0.0e+00  0 84  0  0  0   0 84  0  0  0  8674
MatILUFactorSym        1 1.0 7.8070e-03 2.4 0.00e+00 0.0 0.0e+00 0.0e+00 3.0e+00  1  0  0  0  5   1  0  0  0  5     0
MatAssemblyBegin       2 1.0 1.5763e-0274.4 0.00e+00 0.0 8.4e+01 4.2e+03 4.0e+00  2  0 16 64  7   2  0 16 64  7     0
MatAssemblyEnd         2 1.0 5.6047e-03 1.0 0.00e+00 0.0 1.0e+02 2.2e+02 8.0e+00  1  0 19  4 14   1  0 19  4 15     0
MatGetRowIJ            1 1.0 1.5020e-05 2.5 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatGetOrdering         1 1.0 1.7309e-04 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 4.0e+00  0  0  0  0  7   0  0  0  0  7     0
MatZeroEntries         3 1.0 7.7724e-05 1.5 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPSetUp               2 1.0 1.4210e-04 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPSolve               1 1.0 2.1334e-02 1.0 3.06e+06 3.0 1.5e+02 8.6e+02 1.8e+01  4100 29 24 32   4100 29 24 33  1221
PCSetUp                2 1.0 1.1165e-02 2.2 2.64e+06 3.3 0.0e+00 0.0e+00 9.0e+00  2 84  0  0 16   2 84  0  0 16  1955
PCSetUpOnBlocks        1 1.0 1.0671e-02 2.4 2.64e+06 3.3 0.0e+00 0.0e+00 7.0e+00  1 84  0  0 12   1 84  0  0 13  2046
PCApply                4 1.0 4.3583e-04 1.4 3.36e+05 2.0 0.0e+00 0.0e+00 0.0e+00  0 13  0  0  0   0 13  0  0  0  7510
------------------------------------------------------------------------------------------------------------------------

Memory usage is given in bytes:

Object Type          Creations   Destructions     Memory  Descendants' Mem.
Reports information only for process 0.

--- Event Stage 0: Main Stage

              Vector    12             12        42816     0
      Vector Scatter     2              2         2072     0
           Index Set     9              9        12008     0
   IS L to G Mapping     1              1          564     0
              Matrix     4              4       618248     0
       Krylov Solver     2              2         2368     0
      Preconditioner     2              2         1784     0
              Viewer     1              0            0     0
========================================================================================================================
Average time to get PetscTime(): 9.53674e-08
Average time for MPI_Barrier(): 4.00543e-06
Average time for zero size MPI_Send(): 1.31726e-05
#PETSc Option Table entries:
-ksp_right_pc
-ksp_type cg
-log_summary
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
| Time:           Thu Jan 31 22:15:13 2013                                                                             |
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
| libMesh Performance: Alive time=0.716958, Active time=0.547738                                                 |
 ----------------------------------------------------------------------------------------------------------------
| Event                              nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                              w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|----------------------------------------------------------------------------------------------------------------|
|                                                                                                                |
|                                                                                                                |
| DofMap                                                                                                         |
|   SCALAR_dof_indices()             210       0.0015      0.000007    0.0015      0.000007    0.27     0.27     |
|   add_neighbors_to_send_list()     1         0.0125      0.012467    0.0252      0.025231    2.28     4.61     |
|   build_constraint_matrix()        42        0.0004      0.000010    0.0004      0.000010    0.07     0.07     |
|   build_sparsity()                 1         0.0107      0.010685    0.0266      0.026630    1.95     4.86     |
|   cnstrn_elem_mat_vec()            42        0.0001      0.000002    0.0001      0.000002    0.02     0.02     |
|   create_dof_constraints()         1         0.0191      0.019139    0.1657      0.165707    3.49     30.25    |
|   distribute_dofs()                1         0.0327      0.032680    0.0812      0.081181    5.97     14.82    |
|   dof_indices()                    1378      0.2044      0.000148    0.2060      0.000149    37.32    37.60    |
|   prepare_send_list()              1         0.0002      0.000203    0.0002      0.000203    0.04     0.04     |
|   reinit()                         1         0.0468      0.046847    0.0468      0.046847    8.55     8.55     |
|                                                                                                                |
| EquationSystems                                                                                                |
|   build_solution_vector()          1         0.0019      0.001936    0.0153      0.015340    0.35     2.80     |
|                                                                                                                |
| ExodusII_IO                                                                                                    |
|   write_nodal_data()               1         0.0079      0.007866    0.0079      0.007866    1.44     1.44     |
|                                                                                                                |
| FE                                                                                                             |
|   compute_shape_functions()        71        0.0014      0.000019    0.0014      0.000019    0.25     0.25     |
|   init_shape_functions()           30        0.0003      0.000011    0.0003      0.000011    0.06     0.06     |
|   inverse_map()                    87        0.0015      0.000017    0.0015      0.000017    0.27     0.27     |
|                                                                                                                |
| FEMap                                                                                                          |
|   compute_affine_map()             71        0.0012      0.000017    0.0012      0.000017    0.22     0.22     |
|   compute_face_map()               29        0.0011      0.000039    0.0026      0.000090    0.21     0.48     |
|   init_face_shape_functions()      21        0.0002      0.000011    0.0002      0.000011    0.04     0.04     |
|   init_reference_to_physical_map() 30        0.0013      0.000044    0.0013      0.000044    0.24     0.24     |
|                                                                                                                |
| Mesh                                                                                                           |
|   find_neighbors()                 1         0.0106      0.010596    0.0115      0.011517    1.93     2.10     |
|   renumber_nodes_and_elem()        2         0.0011      0.000561    0.0011      0.000561    0.20     0.20     |
|                                                                                                                |
| MeshCommunication                                                                                              |
|   compute_hilbert_indices()        2         0.0085      0.004240    0.0085      0.004240    1.55     1.55     |
|   find_global_indices()            2         0.0036      0.001781    0.0168      0.008404    0.65     3.07     |
|   parallel_sort()                  2         0.0029      0.001459    0.0035      0.001733    0.53     0.63     |
|                                                                                                                |
| MeshOutput                                                                                                     |
|   write_equation_systems()         1         0.0001      0.000127    0.0235      0.023461    0.02     4.28     |
|                                                                                                                |
| MeshTools::Generation                                                                                          |
|   build_cube()                     1         0.0056      0.005557    0.0056      0.005557    1.01     1.01     |
|                                                                                                                |
| MetisPartitioner                                                                                               |
|   partition()                      1         0.0360      0.035989    0.0436      0.043586    6.57     7.96     |
|                                                                                                                |
| Parallel                                                                                                       |
|   allgather()                      9         0.0006      0.000062    0.0006      0.000069    0.10     0.11     |
|   max(bool)                        1         0.0000      0.000008    0.0000      0.000008    0.00     0.00     |
|   max(scalar)                      105       0.0007      0.000007    0.0007      0.000007    0.14     0.14     |
|   max(vector)                      24        0.0003      0.000013    0.0008      0.000032    0.06     0.14     |
|   min(bool)                        121       0.0007      0.000006    0.0007      0.000006    0.14     0.14     |
|   min(scalar)                      99        0.0122      0.000123    0.0122      0.000123    2.23     2.23     |
|   min(vector)                      24        0.0004      0.000017    0.0010      0.000042    0.08     0.19     |
|   probe()                          132       0.0014      0.000010    0.0014      0.000010    0.25     0.25     |
|   receive()                        132       0.0009      0.000007    0.0023      0.000017    0.16     0.42     |
|   send()                           132       0.0004      0.000003    0.0004      0.000003    0.08     0.08     |
|   send_receive()                   136       0.0012      0.000009    0.0043      0.000032    0.22     0.79     |
|   sum()                            20        0.0008      0.000042    0.0012      0.000060    0.15     0.22     |
|                                                                                                                |
| Parallel::Request                                                                                              |
|   wait()                           132       0.0003      0.000002    0.0003      0.000002    0.05     0.05     |
|                                                                                                                |
| Partitioner                                                                                                    |
|   set_node_processor_ids()         1         0.0021      0.002110    0.0032      0.003235    0.39     0.59     |
|   set_parent_processor_ids()       1         0.0010      0.000984    0.0010      0.000984    0.18     0.18     |
|                                                                                                                |
| PetscLinearSolver                                                                                              |
|   solve()                          1         0.0405      0.040456    0.0405      0.040456    7.39     7.39     |
|                                                                                                                |
| System                                                                                                         |
|   assemble()                       1         0.0706      0.070553    0.1000      0.100046    12.88    18.27    |
 ----------------------------------------------------------------------------------------------------------------
| Totals:                            3102      0.5477                                          100.00            |
 ----------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example systems_of_equations_ex5:
*  mpirun -np 12 example-devel -ksp_type cg -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
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
