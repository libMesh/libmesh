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
        
</pre>
</div>
<div class = "comment">
Create a 2D mesh distributed across the default MPI communicator.
</div>

<div class ="fragment">
<pre>
          Mesh mesh(init.comm(), dim);
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
          system.add_variable("lambda", FIRST, SCALAR);
        
        
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
  
    Mesh mesh(init.comm(), dim);
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
  
    system.add_variable(<B><FONT COLOR="#BC8F8F">&quot;lambda&quot;</FONT></B>, FIRST, SCALAR);
  
  
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
make[4]: Entering directory `/net/spark/workspace/roystgnr/libmesh/git/devel/examples/systems_of_equations/systems_of_equations_ex5'
***************************************************************
* Running Example systems_of_equations_ex5:
*  mpirun -np 4 example-devel -ksp_type cg -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc
***************************************************************
 
 Mesh Information:
  mesh_dimension()=2
  spatial_dimension()=3
  n_nodes()=2121
    n_local_nodes()=547
  n_elem()=500
    n_local_elem()=125
    n_active_elem()=500
  n_subdomains()=1
  n_partitions()=4
  n_processors()=4
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
    n_local_dofs()=1094
    n_constrained_dofs()=42
    n_local_constrained_dofs()=0
    n_vectors()=1
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 30.3929
      Average Off-Processor Bandwidth <= 2.50342
      Maximum  On-Processor Bandwidth <= 1007
      Maximum Off-Processor Bandwidth <= 3236
    DofMap Constraints
      Number of DoF Constraints = 42
      Average DoF Constraint Length= 0
      Number of Node Constraints = 0


 -------------------------------------------------------------------------------------------------------------------
| Processor id:   0                                                                                                 |
| Num Processors: 4                                                                                                 |
| Time:           Fri Apr 19 11:54:15 2013                                                                          |
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
| libMesh Performance: Alive time=0.324504, Active time=0.312767                                                 |
 ----------------------------------------------------------------------------------------------------------------
| Event                              nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                              w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|----------------------------------------------------------------------------------------------------------------|
|                                                                                                                |
|                                                                                                                |
| DofMap                                                                                                         |
|   SCALAR_dof_indices()             522       0.0005      0.000001    0.0005      0.000001    0.15     0.15     |
|   add_neighbors_to_send_list()     1         0.0015      0.001480    0.0018      0.001790    0.47     0.57     |
|   build_constraint_matrix()        125       0.0001      0.000001    0.0001      0.000001    0.04     0.04     |
|   build_sparsity()                 1         0.0015      0.001458    0.0062      0.006175    0.47     1.97     |
|   cnstrn_elem_mat_vec()            125       0.0001      0.000001    0.0001      0.000001    0.03     0.03     |
|   create_dof_constraints()         1         0.0019      0.001949    0.0090      0.009011    0.62     2.88     |
|   distribute_dofs()                1         0.0022      0.002199    0.0097      0.009672    0.70     3.09     |
|   dof_indices()                    2022      0.0134      0.000007    0.0139      0.000007    4.29     4.45     |
|   prepare_send_list()              1         0.0000      0.000011    0.0000      0.000011    0.00     0.00     |
|   reinit()                         1         0.0028      0.002754    0.0028      0.002754    0.88     0.88     |
|                                                                                                                |
| EquationSystems                                                                                                |
|   build_solution_vector()          1         0.0007      0.000687    0.0050      0.004991    0.22     1.60     |
|                                                                                                                |
| ExodusII_IO                                                                                                    |
|   write_nodal_data()               1         0.0988      0.098761    0.0988      0.098761    31.58    31.58    |
|                                                                                                                |
| FE                                                                                                             |
|   compute_shape_functions()        180       0.0005      0.000003    0.0005      0.000003    0.16     0.16     |
|   init_shape_functions()           56        0.0001      0.000002    0.0001      0.000002    0.03     0.03     |
|   inverse_map()                    165       0.0006      0.000004    0.0006      0.000004    0.20     0.20     |
|                                                                                                                |
| FEMap                                                                                                          |
|   compute_affine_map()             180       0.0005      0.000003    0.0005      0.000003    0.17     0.17     |
|   compute_face_map()               55        0.0004      0.000007    0.0010      0.000019    0.12     0.33     |
|   init_face_shape_functions()      21        0.0001      0.000002    0.0001      0.000002    0.02     0.02     |
|   init_reference_to_physical_map() 56        0.0004      0.000008    0.0004      0.000008    0.14     0.14     |
|                                                                                                                |
| Mesh                                                                                                           |
|   find_neighbors()                 1         0.0008      0.000813    0.0015      0.001512    0.26     0.48     |
|   renumber_nodes_and_elem()        2         0.0002      0.000094    0.0002      0.000094    0.06     0.06     |
|                                                                                                                |
| MeshCommunication                                                                                              |
|   compute_hilbert_indices()        2         0.0026      0.001311    0.0026      0.001311    0.84     0.84     |
|   find_global_indices()            2         0.0004      0.000203    0.0058      0.002903    0.13     1.86     |
|   parallel_sort()                  2         0.0003      0.000166    0.0024      0.001203    0.11     0.77     |
|                                                                                                                |
| MeshOutput                                                                                                     |
|   write_equation_systems()         1         0.0001      0.000059    0.1039      0.103879    0.02     33.21    |
|                                                                                                                |
| MeshTools::Generation                                                                                          |
|   build_cube()                     1         0.0006      0.000613    0.0006      0.000613    0.20     0.20     |
|                                                                                                                |
| MetisPartitioner                                                                                               |
|   partition()                      1         0.0031      0.003104    0.0063      0.006319    0.99     2.02     |
|                                                                                                                |
| Parallel                                                                                                       |
|   allgather()                      9         0.0031      0.000342    0.0031      0.000349    0.98     1.00     |
|   max(bool)                        1         0.0000      0.000004    0.0000      0.000004    0.00     0.00     |
|   max(scalar)                      105       0.0009      0.000008    0.0009      0.000008    0.28     0.28     |
|   max(vector)                      24        0.0002      0.000009    0.0008      0.000032    0.07     0.25     |
|   min(bool)                        121       0.0009      0.000007    0.0009      0.000007    0.27     0.27     |
|   min(scalar)                      99        0.0204      0.000206    0.0204      0.000206    6.51     6.51     |
|   min(vector)                      24        0.0003      0.000011    0.0014      0.000057    0.09     0.44     |
|   probe()                          36        0.0013      0.000037    0.0013      0.000037    0.43     0.43     |
|   receive()                        36        0.0001      0.000003    0.0015      0.000040    0.03     0.46     |
|   send()                           36        0.0001      0.000002    0.0001      0.000002    0.03     0.03     |
|   send_receive()                   40        0.0002      0.000005    0.0018      0.000044    0.06     0.56     |
|   sum()                            20        0.0023      0.000117    0.0044      0.000219    0.75     1.40     |
|                                                                                                                |
| Parallel::Request                                                                                              |
|   wait()                           36        0.0000      0.000001    0.0000      0.000001    0.01     0.01     |
|                                                                                                                |
| Partitioner                                                                                                    |
|   set_node_processor_ids()         1         0.0005      0.000541    0.0029      0.002887    0.17     0.92     |
|   set_parent_processor_ids()       1         0.0001      0.000085    0.0001      0.000085    0.03     0.03     |
|                                                                                                                |
| PetscLinearSolver                                                                                              |
|   solve()                          1         0.1161      0.116108    0.1161      0.116108    37.12    37.12    |
|                                                                                                                |
| System                                                                                                         |
|   assemble()                       1         0.0321      0.032097    0.0381      0.038065    10.26    12.17    |
 ----------------------------------------------------------------------------------------------------------------
| Totals:                            4118      0.3128                                          100.00            |
 ----------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example systems_of_equations_ex5:
*  mpirun -np 4 example-devel -ksp_type cg -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc
***************************************************************
make[4]: Leaving directory `/net/spark/workspace/roystgnr/libmesh/git/devel/examples/systems_of_equations/systems_of_equations_ex5'
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
