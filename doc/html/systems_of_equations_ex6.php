<?php $root=""; ?>
<?php require($root."navigation.php"); ?>
<html>
<head>
  <?php load_style($root); ?>
</head>
 
<body>
 
<?php make_navigation("systems_of_equations_ex6",$root)?>
 
<div class="content">
<a name="comments"></a> 
<br><br><br> <h1> The source file systems_of_equations_ex6.C with comments: </h1> 
<div class = "comment">
<h1> Systems Example 6 - 3D Linear Elastic Cantilever </h1>
By David Knezevic

<br><br>This is a 3D version of systems_of_equations_ex4.
In this case we also compute and plot stresses.


<br><br>

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
        #include "libmesh/gnuplot_io.h"
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
        
        #define x_scaling 1.3
        #define x_load 0.0
        #define y_load 0.0
        #define z_load -1.0
        
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
Matrix and right-hand side assembly
</div>

<div class ="fragment">
<pre>
        void assemble_elasticity(EquationSystems& es,
                                 const std::string& system_name);
                                 
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
The Kronecker delta function, used in eval_elasticity_tensor
</div>

<div class ="fragment">
<pre>
        Real kronecker_delta(unsigned int i,
                             unsigned int j);
        
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
          const unsigned int dim = 3;
        
</pre>
</div>
<div class = "comment">
Make sure libMesh was compiled for 3D
</div>

<div class ="fragment">
<pre>
          libmesh_example_assert(dim == LIBMESH_DIM, "3D support");
        
          Mesh mesh(dim);
          MeshTools::Generation::build_cube (mesh,
                                             40,
                                             8,
                                             4,
                                             0., 1.*x_scaling,
                                             0., 0.3,
                                             0., 0.1,
                                             HEX8);
        
          
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
Add three displacement variables, u and v, to the system
</div>

<div class ="fragment">
<pre>
          unsigned int u_var = system.add_variable("u", FIRST, LAGRANGE);
          unsigned int v_var = system.add_variable("v", FIRST, LAGRANGE);
          unsigned int w_var = system.add_variable("w", FIRST, LAGRANGE);
        
          system.attach_assemble_function (assemble_elasticity);
        
        
          std::set&lt;boundary_id_type&gt; boundary_ids;
          boundary_ids.insert(BOUNDARY_ID_MIN_X);
        
</pre>
</div>
<div class = "comment">
Create a vector storing the variable numbers which the BC applies to
</div>

<div class ="fragment">
<pre>
          std::vector&lt;unsigned int&gt; variables;
          variables.push_back(u_var);
          variables.push_back(v_var);
          variables.push_back(w_var);
          
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
Post-process the solution to compute the stresses
</div>

<div class ="fragment">
<pre>
          compute_stresses(equation_systems);
        
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
        
          const unsigned int n_components = 3;
          const unsigned int u_var = system.variable_number ("u");
          const unsigned int v_var = system.variable_number ("v");
          const unsigned int w_var = system.variable_number ("w");
        
          const DofMap& dof_map = system.get_dof_map();
          FEType fe_type = dof_map.variable_type(u_var);
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
            Kuu(Ke), Kuv(Ke), Kuw(Ke),
            Kvu(Ke), Kvv(Ke), Kvw(Ke),
            Kwu(Ke), Kwv(Ke), Kww(Ke);
        
          DenseSubVector&lt;Number&gt;
            Fu(Fe),
            Fv(Fe),
            Fw(Fe);
        
          std::vector&lt;dof_id_type&gt; dof_indices;
          std::vector&lt;dof_id_type&gt; dof_indices_u;
          std::vector&lt;dof_id_type&gt; dof_indices_v;
          std::vector&lt;dof_id_type&gt; dof_indices_w;
        
          MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
          const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
        
          for ( ; el != end_el; ++el)
            {
              const Elem* elem = *el;
        
              dof_map.dof_indices (elem, dof_indices);
              dof_map.dof_indices (elem, dof_indices_u, u_var);
              dof_map.dof_indices (elem, dof_indices_v, v_var);
              dof_map.dof_indices (elem, dof_indices_w, w_var);
        
              const unsigned int n_dofs   = dof_indices.size();
              const unsigned int n_u_dofs = dof_indices_u.size(); 
              const unsigned int n_v_dofs = dof_indices_v.size();
              const unsigned int n_w_dofs = dof_indices_w.size();
        
              fe-&gt;reinit (elem);
        
              Ke.resize (n_dofs, n_dofs);
              Fe.resize (n_dofs);
        
              Kuu.reposition (u_var*n_u_dofs, u_var*n_u_dofs, n_u_dofs, n_u_dofs);
              Kuv.reposition (u_var*n_u_dofs, v_var*n_u_dofs, n_u_dofs, n_v_dofs);
              Kuw.reposition (u_var*n_u_dofs, w_var*n_u_dofs, n_u_dofs, n_w_dofs);
              
              Kvu.reposition (v_var*n_u_dofs, u_var*n_u_dofs, n_v_dofs, n_u_dofs);
              Kvv.reposition (v_var*n_u_dofs, v_var*n_u_dofs, n_v_dofs, n_v_dofs);
              Kvw.reposition (v_var*n_u_dofs, w_var*n_u_dofs, n_v_dofs, n_w_dofs);
        
              Kwu.reposition (w_var*n_u_dofs, u_var*n_u_dofs, n_w_dofs, n_u_dofs);
              Kwv.reposition (w_var*n_u_dofs, v_var*n_u_dofs, n_w_dofs, n_v_dofs);
              Kww.reposition (w_var*n_u_dofs, w_var*n_u_dofs, n_w_dofs, n_w_dofs);
              
              Fu.reposition (u_var*n_u_dofs, n_u_dofs);
              Fv.reposition (v_var*n_u_dofs, n_v_dofs);
              Fw.reposition (w_var*n_u_dofs, n_w_dofs);
        
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
                      unsigned int C_i=0, C_k=0;
        
                      for(unsigned int C_j=0; C_j&lt;n_components; C_j++)
                        for(unsigned int C_l=0; C_l&lt;n_components; C_l++)
                        {
                          Kuu(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i,C_j,C_k,C_l) * dphi[i][qp](C_j)*dphi[j][qp](C_l));
                        }
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
                      unsigned int C_i=0, C_k=1;
        
                      for(unsigned int C_j=0; C_j&lt;n_components; C_j++)
                        for(unsigned int C_l=0; C_l&lt;n_components; C_l++)
                        {
                          Kuv(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i,C_j,C_k,C_l) * dphi[i][qp](C_j)*dphi[j][qp](C_l));
                        }
                    }
        
                  for (unsigned int i=0; i&lt;n_u_dofs; i++)
                    for (unsigned int j=0; j&lt;n_w_dofs; j++)
                    {
</pre>
</div>
<div class = "comment">
Tensor indices
</div>

<div class ="fragment">
<pre>
                      unsigned int C_i=0, C_k=2;
        
                      for(unsigned int C_j=0; C_j&lt;n_components; C_j++)
                        for(unsigned int C_l=0; C_l&lt;n_components; C_l++)
                        {
                          Kuw(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i,C_j,C_k,C_l) * dphi[i][qp](C_j)*dphi[j][qp](C_l));
                        }
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
                      unsigned int C_i = 1, C_k = 0;
        
                      for(unsigned int C_j=0; C_j&lt;n_components; C_j++)
                        for(unsigned int C_l=0; C_l&lt;n_components; C_l++)
                        {
                          Kvu(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i,C_j,C_k,C_l) * dphi[i][qp](C_j)*dphi[j][qp](C_l));
                        }
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
                      unsigned int C_i = 1, C_k = 1;
        
                      for(unsigned int C_j=0; C_j&lt;n_components; C_j++)
                        for(unsigned int C_l=0; C_l&lt;n_components; C_l++)
                        {
                          Kvv(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i,C_j,C_k,C_l) * dphi[i][qp](C_j)*dphi[j][qp](C_l));
                        }
                    }
                    
                  for (unsigned int i=0; i&lt;n_v_dofs; i++)
                    for (unsigned int j=0; j&lt;n_w_dofs; j++)
                    {
</pre>
</div>
<div class = "comment">
Tensor indices
</div>

<div class ="fragment">
<pre>
                      unsigned int C_i = 1, C_k = 2;
        
                      for(unsigned int C_j=0; C_j&lt;n_components; C_j++)
                        for(unsigned int C_l=0; C_l&lt;n_components; C_l++)
                        {
                          Kvw(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i,C_j,C_k,C_l) * dphi[i][qp](C_j)*dphi[j][qp](C_l));
                        }
                    }
                    
                  for (unsigned int i=0; i&lt;n_w_dofs; i++)
                    for (unsigned int j=0; j&lt;n_u_dofs; j++)
                    {
</pre>
</div>
<div class = "comment">
Tensor indices
</div>

<div class ="fragment">
<pre>
                      unsigned int C_i = 2, C_k = 0;
        
                      for(unsigned int C_j=0; C_j&lt;n_components; C_j++)
                        for(unsigned int C_l=0; C_l&lt;n_components; C_l++)
                        {
                          Kwu(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i,C_j,C_k,C_l) * dphi[i][qp](C_j)*dphi[j][qp](C_l));
                        }
                    }
                    
                  for (unsigned int i=0; i&lt;n_w_dofs; i++)
                    for (unsigned int j=0; j&lt;n_v_dofs; j++)
                    {
</pre>
</div>
<div class = "comment">
Tensor indices
</div>

<div class ="fragment">
<pre>
                      unsigned int C_i = 2, C_k = 1;
        
                      for(unsigned int C_j=0; C_j&lt;n_components; C_j++)
                        for(unsigned int C_l=0; C_l&lt;n_components; C_l++)
                        {
                          Kwv(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i,C_j,C_k,C_l) * dphi[i][qp](C_j)*dphi[j][qp](C_l));
                        }
                    }
                    
                  for (unsigned int i=0; i&lt;n_w_dofs; i++)
                    for (unsigned int j=0; j&lt;n_w_dofs; j++)
                    {
</pre>
</div>
<div class = "comment">
Tensor indices
</div>

<div class ="fragment">
<pre>
                      unsigned int C_i = 2, C_k = 2;
        
                      for(unsigned int C_j=0; C_j&lt;n_components; C_j++)
                        for(unsigned int C_l=0; C_l&lt;n_components; C_l++)
                        {
                          Kww(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i,C_j,C_k,C_l) * dphi[i][qp](C_j)*dphi[j][qp](C_l));
                        }
                    }
              }
        
              {
                for (unsigned int side=0; side&lt;elem-&gt;n_sides(); side++)
                  if (elem-&gt;neighbor(side) == NULL)
                    {
                      const std::vector&lt;std::vector&lt;Real&gt; &gt;&  phi_face = fe_face-&gt;get_phi();
                      const std::vector&lt;Real&gt;& JxW_face = fe_face-&gt;get_JxW();
        
                      fe_face-&gt;reinit(elem, side);
        
                      for (unsigned int qp=0; qp&lt;qface.n_points(); qp++)
                      {
</pre>
</div>
<div class = "comment">
Apply a traction
</div>

<div class ="fragment">
<pre>
                        if( mesh.boundary_info-&gt;has_boundary_id(elem, side, BOUNDARY_ID_MAX_X) )
                        {
                          for (unsigned int i=0; i&lt;n_v_dofs; i++)
                          {
                            Fu(i) += JxW_face[qp] * x_load * phi_face[i][qp];
                          }
                          for (unsigned int i=0; i&lt;n_v_dofs; i++)
                          {
                            Fv(i) += JxW_face[qp] * y_load * phi_face[i][qp];
                          }
                          for (unsigned int i=0; i&lt;n_v_dofs; i++)
                          {
                            Fw(i) += JxW_face[qp] * z_load * phi_face[i][qp];
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
        
        void compute_stresses(EquationSystems& es)
        {
          const MeshBase& mesh = es.get_mesh();
        
          const unsigned int dim = mesh.mesh_dimension();
        
          LinearImplicitSystem& system = es.get_system&lt;LinearImplicitSystem&gt;("Elasticity");
        
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
                      elem_sigma(C_i,C_j) += JxW[qp]*( eval_elasticity_tensor(C_i,C_j,C_k,C_l) * displacement_gradient(C_l) );
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
        }
        
        Real kronecker_delta(unsigned int i,
                             unsigned int j)
        {
          return i == j ? 1. : 0.;
        }
        
        Real eval_elasticity_tensor(unsigned int i,
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
</pre>
</div>

<a name="nocomments"></a> 
<br><br><br> <h1> The source file systems_of_equations_ex6.C without comments: </h1> 
<pre> 
  
  
  #include &lt;iostream&gt;
  #include &lt;algorithm&gt;
  #include &lt;math.h&gt;
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/libmesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/mesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/mesh_generation.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/exodusII_io.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/gnuplot_io.h&quot;</FONT></B>
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
  
  #define x_scaling 1.3
  #define x_load 0.0
  #define y_load 0.0
  #define z_load -1.0
  
  #define BOUNDARY_ID_MIN_Z 0
  #define BOUNDARY_ID_MIN_Y 1
  #define BOUNDARY_ID_MAX_X 2
  #define BOUNDARY_ID_MAX_Y 3
  #define BOUNDARY_ID_MIN_X 4
  #define BOUNDARY_ID_MAX_Z 5
  
  using namespace libMesh;
  
  <B><FONT COLOR="#228B22">void</FONT></B> assemble_elasticity(EquationSystems&amp; es,
                           <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp; system_name);
                           
  <B><FONT COLOR="#228B22">void</FONT></B> compute_stresses(EquationSystems&amp; es);
  
  Real kronecker_delta(<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i,
                       <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j);
  
  Real eval_elasticity_tensor(<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i,
                              <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j,
                              <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> k,
                              <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> l);
  
  <B><FONT COLOR="#228B22">int</FONT></B> main (<B><FONT COLOR="#228B22">int</FONT></B> argc, <B><FONT COLOR="#228B22">char</FONT></B>** argv)
  {
    LibMeshInit init (argc, argv);
  
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> dim = 3;
  
    libmesh_example_assert(dim == LIBMESH_DIM, <B><FONT COLOR="#BC8F8F">&quot;3D support&quot;</FONT></B>);
  
    Mesh mesh(dim);
    <B><FONT COLOR="#5F9EA0">MeshTools</FONT></B>::Generation::build_cube (mesh,
                                       40,
                                       8,
                                       4,
                                       0., 1.*x_scaling,
                                       0., 0.3,
                                       0., 0.1,
                                       HEX8);
  
    
    mesh.print_info();
    
    
    EquationSystems equation_systems (mesh);
    
    LinearImplicitSystem&amp; system =
      equation_systems.add_system&lt;LinearImplicitSystem&gt; (<B><FONT COLOR="#BC8F8F">&quot;Elasticity&quot;</FONT></B>);
    
    <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> u_var = system.add_variable(<B><FONT COLOR="#BC8F8F">&quot;u&quot;</FONT></B>, FIRST, LAGRANGE);
    <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> v_var = system.add_variable(<B><FONT COLOR="#BC8F8F">&quot;v&quot;</FONT></B>, FIRST, LAGRANGE);
    <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> w_var = system.add_variable(<B><FONT COLOR="#BC8F8F">&quot;w&quot;</FONT></B>, FIRST, LAGRANGE);
  
    system.attach_assemble_function (assemble_elasticity);
  
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::set&lt;boundary_id_type&gt; boundary_ids;
    boundary_ids.insert(BOUNDARY_ID_MIN_X);
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B>&gt; variables;
    variables.push_back(u_var);
    variables.push_back(v_var);
    variables.push_back(w_var);
    
    ZeroFunction&lt;&gt; zf;
    
    DirichletBoundary dirichlet_bc(boundary_ids,
                                   variables,
                                   &amp;zf);
  
    system.get_dof_map().add_dirichlet_boundary(dirichlet_bc);
    
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
    
    equation_systems.init();
  
    equation_systems.print_info();
  
    system.solve();
    
    compute_stresses(equation_systems);
  
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
  
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_components = 3;
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> u_var = system.variable_number (<B><FONT COLOR="#BC8F8F">&quot;u&quot;</FONT></B>);
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> v_var = system.variable_number (<B><FONT COLOR="#BC8F8F">&quot;v&quot;</FONT></B>);
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> w_var = system.variable_number (<B><FONT COLOR="#BC8F8F">&quot;w&quot;</FONT></B>);
  
    <B><FONT COLOR="#228B22">const</FONT></B> DofMap&amp; dof_map = system.get_dof_map();
    FEType fe_type = dof_map.variable_type(u_var);
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
      Kuu(Ke), Kuv(Ke), Kuw(Ke),
      Kvu(Ke), Kvv(Ke), Kvw(Ke),
      Kwu(Ke), Kwv(Ke), Kww(Ke);
  
    DenseSubVector&lt;Number&gt;
      Fu(Fe),
      Fv(Fe),
      Fw(Fe);
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;dof_id_type&gt; dof_indices;
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;dof_id_type&gt; dof_indices_u;
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;dof_id_type&gt; dof_indices_v;
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;dof_id_type&gt; dof_indices_w;
  
    <B><FONT COLOR="#5F9EA0">MeshBase</FONT></B>::const_element_iterator       el     = mesh.active_local_elements_begin();
    <B><FONT COLOR="#228B22">const</FONT></B> MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
  
    <B><FONT COLOR="#A020F0">for</FONT></B> ( ; el != end_el; ++el)
      {
        <B><FONT COLOR="#228B22">const</FONT></B> Elem* elem = *el;
  
        dof_map.dof_indices (elem, dof_indices);
        dof_map.dof_indices (elem, dof_indices_u, u_var);
        dof_map.dof_indices (elem, dof_indices_v, v_var);
        dof_map.dof_indices (elem, dof_indices_w, w_var);
  
        <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_dofs   = dof_indices.size();
        <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_u_dofs = dof_indices_u.size(); 
        <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_v_dofs = dof_indices_v.size();
        <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_w_dofs = dof_indices_w.size();
  
        fe-&gt;reinit (elem);
  
        Ke.resize (n_dofs, n_dofs);
        Fe.resize (n_dofs);
  
        Kuu.reposition (u_var*n_u_dofs, u_var*n_u_dofs, n_u_dofs, n_u_dofs);
        Kuv.reposition (u_var*n_u_dofs, v_var*n_u_dofs, n_u_dofs, n_v_dofs);
        Kuw.reposition (u_var*n_u_dofs, w_var*n_u_dofs, n_u_dofs, n_w_dofs);
        
        Kvu.reposition (v_var*n_u_dofs, u_var*n_u_dofs, n_v_dofs, n_u_dofs);
        Kvv.reposition (v_var*n_u_dofs, v_var*n_u_dofs, n_v_dofs, n_v_dofs);
        Kvw.reposition (v_var*n_u_dofs, w_var*n_u_dofs, n_v_dofs, n_w_dofs);
  
        Kwu.reposition (w_var*n_u_dofs, u_var*n_u_dofs, n_w_dofs, n_u_dofs);
        Kwv.reposition (w_var*n_u_dofs, v_var*n_u_dofs, n_w_dofs, n_v_dofs);
        Kww.reposition (w_var*n_u_dofs, w_var*n_u_dofs, n_w_dofs, n_w_dofs);
        
        Fu.reposition (u_var*n_u_dofs, n_u_dofs);
        Fv.reposition (v_var*n_u_dofs, n_v_dofs);
        Fw.reposition (w_var*n_u_dofs, n_w_dofs);
  
        <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> qp=0; qp&lt;qrule.n_points(); qp++)
        {
            <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i&lt;n_u_dofs; i++)
              <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j=0; j&lt;n_u_dofs; j++)
              {
                <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> C_i=0, C_k=0;
  
                <B><FONT COLOR="#A020F0">for</FONT></B>(<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> C_j=0; C_j&lt;n_components; C_j++)
                  <B><FONT COLOR="#A020F0">for</FONT></B>(<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> C_l=0; C_l&lt;n_components; C_l++)
                  {
                    Kuu(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i,C_j,C_k,C_l) * dphi[i][qp](C_j)*dphi[j][qp](C_l));
                  }
              }
  
            <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i&lt;n_u_dofs; i++)
              <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j=0; j&lt;n_v_dofs; j++)
              {
                <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> C_i=0, C_k=1;
  
                <B><FONT COLOR="#A020F0">for</FONT></B>(<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> C_j=0; C_j&lt;n_components; C_j++)
                  <B><FONT COLOR="#A020F0">for</FONT></B>(<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> C_l=0; C_l&lt;n_components; C_l++)
                  {
                    Kuv(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i,C_j,C_k,C_l) * dphi[i][qp](C_j)*dphi[j][qp](C_l));
                  }
              }
  
            <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i&lt;n_u_dofs; i++)
              <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j=0; j&lt;n_w_dofs; j++)
              {
                <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> C_i=0, C_k=2;
  
                <B><FONT COLOR="#A020F0">for</FONT></B>(<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> C_j=0; C_j&lt;n_components; C_j++)
                  <B><FONT COLOR="#A020F0">for</FONT></B>(<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> C_l=0; C_l&lt;n_components; C_l++)
                  {
                    Kuw(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i,C_j,C_k,C_l) * dphi[i][qp](C_j)*dphi[j][qp](C_l));
                  }
              }
  
            <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i&lt;n_v_dofs; i++)
              <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j=0; j&lt;n_u_dofs; j++)
              {
                <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> C_i = 1, C_k = 0;
  
                <B><FONT COLOR="#A020F0">for</FONT></B>(<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> C_j=0; C_j&lt;n_components; C_j++)
                  <B><FONT COLOR="#A020F0">for</FONT></B>(<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> C_l=0; C_l&lt;n_components; C_l++)
                  {
                    Kvu(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i,C_j,C_k,C_l) * dphi[i][qp](C_j)*dphi[j][qp](C_l));
                  }
              }
  
            <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i&lt;n_v_dofs; i++)
              <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j=0; j&lt;n_v_dofs; j++)
              {
                <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> C_i = 1, C_k = 1;
  
                <B><FONT COLOR="#A020F0">for</FONT></B>(<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> C_j=0; C_j&lt;n_components; C_j++)
                  <B><FONT COLOR="#A020F0">for</FONT></B>(<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> C_l=0; C_l&lt;n_components; C_l++)
                  {
                    Kvv(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i,C_j,C_k,C_l) * dphi[i][qp](C_j)*dphi[j][qp](C_l));
                  }
              }
              
            <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i&lt;n_v_dofs; i++)
              <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j=0; j&lt;n_w_dofs; j++)
              {
                <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> C_i = 1, C_k = 2;
  
                <B><FONT COLOR="#A020F0">for</FONT></B>(<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> C_j=0; C_j&lt;n_components; C_j++)
                  <B><FONT COLOR="#A020F0">for</FONT></B>(<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> C_l=0; C_l&lt;n_components; C_l++)
                  {
                    Kvw(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i,C_j,C_k,C_l) * dphi[i][qp](C_j)*dphi[j][qp](C_l));
                  }
              }
              
            <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i&lt;n_w_dofs; i++)
              <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j=0; j&lt;n_u_dofs; j++)
              {
                <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> C_i = 2, C_k = 0;
  
                <B><FONT COLOR="#A020F0">for</FONT></B>(<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> C_j=0; C_j&lt;n_components; C_j++)
                  <B><FONT COLOR="#A020F0">for</FONT></B>(<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> C_l=0; C_l&lt;n_components; C_l++)
                  {
                    Kwu(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i,C_j,C_k,C_l) * dphi[i][qp](C_j)*dphi[j][qp](C_l));
                  }
              }
              
            <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i&lt;n_w_dofs; i++)
              <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j=0; j&lt;n_v_dofs; j++)
              {
                <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> C_i = 2, C_k = 1;
  
                <B><FONT COLOR="#A020F0">for</FONT></B>(<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> C_j=0; C_j&lt;n_components; C_j++)
                  <B><FONT COLOR="#A020F0">for</FONT></B>(<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> C_l=0; C_l&lt;n_components; C_l++)
                  {
                    Kwv(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i,C_j,C_k,C_l) * dphi[i][qp](C_j)*dphi[j][qp](C_l));
                  }
              }
              
            <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i&lt;n_w_dofs; i++)
              <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j=0; j&lt;n_w_dofs; j++)
              {
                <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> C_i = 2, C_k = 2;
  
                <B><FONT COLOR="#A020F0">for</FONT></B>(<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> C_j=0; C_j&lt;n_components; C_j++)
                  <B><FONT COLOR="#A020F0">for</FONT></B>(<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> C_l=0; C_l&lt;n_components; C_l++)
                  {
                    Kww(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i,C_j,C_k,C_l) * dphi[i][qp](C_j)*dphi[j][qp](C_l));
                  }
              }
        }
  
        {
          <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> side=0; side&lt;elem-&gt;n_sides(); side++)
            <B><FONT COLOR="#A020F0">if</FONT></B> (elem-&gt;neighbor(side) == NULL)
              {
                <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;Real&gt; &gt;&amp;  phi_face = fe_face-&gt;get_phi();
                <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Real&gt;&amp; JxW_face = fe_face-&gt;get_JxW();
  
                fe_face-&gt;reinit(elem, side);
  
                <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> qp=0; qp&lt;qface.n_points(); qp++)
                {
                  <B><FONT COLOR="#A020F0">if</FONT></B>( mesh.boundary_info-&gt;has_boundary_id(elem, side, BOUNDARY_ID_MAX_X) )
                  {
                    <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i&lt;n_v_dofs; i++)
                    {
                      Fu(i) += JxW_face[qp] * x_load * phi_face[i][qp];
                    }
                    <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i&lt;n_v_dofs; i++)
                    {
                      Fv(i) += JxW_face[qp] * y_load * phi_face[i][qp];
                    }
                    <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i&lt;n_v_dofs; i++)
                    {
                      Fw(i) += JxW_face[qp] * z_load * phi_face[i][qp];
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
  
  <B><FONT COLOR="#228B22">void</FONT></B> compute_stresses(EquationSystems&amp; es)
  {
    <B><FONT COLOR="#228B22">const</FONT></B> MeshBase&amp; mesh = es.get_mesh();
  
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> dim = mesh.mesh_dimension();
  
    LinearImplicitSystem&amp; system = es.get_system&lt;LinearImplicitSystem&gt;(<B><FONT COLOR="#BC8F8F">&quot;Elasticity&quot;</FONT></B>);
  
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
                elem_sigma(C_i,C_j) += JxW[qp]*( eval_elasticity_tensor(C_i,C_j,C_k,C_l) * displacement_gradient(C_l) );
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
  }
  
  Real kronecker_delta(<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i,
                       <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j)
  {
    <B><FONT COLOR="#A020F0">return</FONT></B> i == j ? 1. : 0.;
  }
  
  Real eval_elasticity_tensor(<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i,
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
</pre> 
<a name="output"></a> 
<br><br><br> <h1> The console output of the program: </h1> 
<pre>
***************************************************************
* Running Example systems_of_equations_ex6:
*  mpirun -np 12 example-devel -ksp_type cg -pc_type bjacobi -sub_pc_type icc -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
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

 EquationSystems
  n_systems()=2
   System #0, "Elasticity"
    Type "LinearImplicit"
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

************************************************************************************************************************
***             WIDEN YOUR WINDOW TO 120 CHARACTERS.  Use 'enscript -r -fCourier9' to print this document            ***
************************************************************************************************************************

---------------------------------------------- PETSc Performance Summary: ----------------------------------------------

/workspace/libmesh/examples/systems_of_equations/systems_of_equations_ex6/.libs/lt-example-devel on a intel-12. named hbar.ices.utexas.edu with 12 processors, by benkirk Thu Jan 31 22:15:43 2013
Using Petsc Release Version 3.3.0, Patch 2, Fri Jul 13 15:42:00 CDT 2012 

                         Max       Max/Min        Avg      Total 
Time (sec):           2.879e+00      1.00003   2.879e+00
Objects:              4.100e+01      1.00000   4.100e+01
Flops:                8.200e+07      2.90808   5.548e+07  6.658e+08
Flops/sec:            2.849e+07      2.90815   1.927e+07  2.313e+08
MPI Messages:         6.825e+02      1.50000   5.688e+02  6.825e+03
MPI Message Lengths:  7.276e+05      1.65309   1.096e+03  7.479e+06
MPI Reductions:       7.130e+02      1.00000

Flop counting convention: 1 flop = 1 real number operation of type (multiply/divide/add/subtract)
                            e.g., VecAXPY() for real vectors of length N --> 2N flops
                            and VecAXPY() for complex vectors of length N --> 8N flops

Summary of Stages:   ----- Time ------  ----- Flops -----  --- Messages ---  -- Message Lengths --  -- Reductions --
                        Avg     %Total     Avg     %Total   counts   %Total     Avg         %Total   counts   %Total 
 0:      Main Stage: 2.8787e+00 100.0%  6.6576e+08 100.0%  6.825e+03 100.0%  1.096e+03      100.0%  7.120e+02  99.9% 

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

VecTDot              432 1.0 9.8159e-03 3.1 5.18e+05 2.0 0.0e+00 0.0e+00 4.3e+02  0  1  0  0 61   0  1  0  0 61   487
VecNorm              218 1.0 5.8289e-02 5.2 2.62e+05 2.0 0.0e+00 0.0e+00 2.2e+02  1  0  0  0 31   1  0  0  0 31    41
VecCopy                3 1.0 9.2983e-06 1.3 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecSet               226 1.0 3.1471e-04 1.3 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecAXPY              432 1.0 8.2636e-04 1.5 5.18e+05 2.0 0.0e+00 0.0e+00 0.0e+00  0  1  0  0  0   0  1  0  0  0  5787
VecAYPX              216 1.0 3.8314e-04 1.6 2.59e+05 2.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0  6226
VecAssemblyBegin       5 1.0 3.9686e-02111.6 0.00e+00 0.0 3.0e+01 2.0e+03 1.5e+01  1  0  0  1  2   1  0  0  1  2     0
VecAssemblyEnd         5 1.0 4.2200e-05 1.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecScatterBegin      219 1.0 1.4424e-03 1.2 0.00e+00 0.0 6.6e+03 8.6e+02 0.0e+00  0  0 96 75  0   0  0 96 75  0     0
VecScatterEnd        219 1.0 3.5646e-0254.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatMult              217 1.0 5.1867e-02 2.4 1.68e+07 2.0 6.5e+03 8.5e+02 0.0e+00  1 23 95 74  0   1 23 95 74  0  2938
MatSolve             218 1.0 5.5080e-02 2.8 5.09e+07 3.1 0.0e+00 0.0e+00 0.0e+00  1 62  0  0  0   1 62  0  0  0  7462
MatLUFactorNum         1 1.0 1.0329e-02 4.6 1.28e+07 4.8 0.0e+00 0.0e+00 0.0e+00  0 13  0  0  0   0 13  0  0  0  8517
MatILUFactorSym        1 1.0 3.5388e-02 3.7 0.00e+00 0.0 0.0e+00 0.0e+00 3.0e+00  1  0  0  0  0   1  0  0  0  0     0
MatAssemblyBegin       2 1.0 8.1562e-01142.0 0.00e+00 0.0 4.5e+01 3.8e+04 4.0e+00 14  0  1 23  1  14  0  1 23  1     0
MatAssemblyEnd         2 1.0 3.4211e-03 1.8 0.00e+00 0.0 6.0e+01 2.1e+02 8.0e+00  0  0  1  0  1   0  0  1  0  1     0
MatGetRowIJ            1 1.0 1.3113e-05 4.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatGetOrdering         1 1.0 1.6499e-04 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 4.0e+00  0  0  0  0  1   0  0  0  0  1     0
MatZeroEntries         3 1.0 2.5344e-04 2.6 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPSetUp               2 1.0 1.5497e-04 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPSolve               1 1.0 1.4923e-01 1.0 8.20e+07 2.9 6.5e+03 8.5e+02 6.6e+02  5100 95 74 92   5100 95 74 93  4461
PCSetUp                2 1.0 4.6556e-02 3.7 1.28e+07 4.8 0.0e+00 0.0e+00 9.0e+00  1 13  0  0  1   1 13  0  0  1  1890
PCSetUpOnBlocks        1 1.0 4.6039e-02 3.8 1.28e+07 4.8 0.0e+00 0.0e+00 7.0e+00  1 13  0  0  1   1 13  0  0  1  1911
PCApply              218 1.0 6.1135e-02 2.4 5.09e+07 3.1 0.0e+00 0.0e+00 0.0e+00  2 62  0  0  0   2 62  0  0  0  6723
------------------------------------------------------------------------------------------------------------------------

Memory usage is given in bytes:

Object Type          Creations   Destructions     Memory  Descendants' Mem.
Reports information only for process 0.

--- Event Stage 0: Main Stage

              Vector    16             16        93136     0
      Vector Scatter     3              3         3108     0
           Index Set    11             11        18760     0
   IS L to G Mapping     2              2         1128     0
              Matrix     4              4      1968640     0
       Krylov Solver     2              2         2368     0
      Preconditioner     2              2         1784     0
              Viewer     1              0            0     0
========================================================================================================================
Average time to get PetscTime(): 1.19209e-07
Average time for MPI_Barrier(): 5.00679e-06
Average time for zero size MPI_Send(): 1.71661e-05
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
| Time:           Thu Jan 31 22:15:43 2013                                                                             |
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
| libMesh Performance: Alive time=2.99182, Active time=2.75701                                                   |
 ----------------------------------------------------------------------------------------------------------------
| Event                              nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                              w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|----------------------------------------------------------------------------------------------------------------|
|                                                                                                                |
|                                                                                                                |
| DofMap                                                                                                         |
|   add_neighbors_to_send_list()     2         0.0602      0.030108    0.1266      0.063312    2.18     4.59     |
|   build_constraint_matrix()        108       0.0013      0.000012    0.0013      0.000012    0.05     0.05     |
|   build_sparsity()                 1         0.0525      0.052483    0.0973      0.097266    1.90     3.53     |
|   cnstrn_elem_mat_vec()            108       0.0002      0.000002    0.0002      0.000002    0.01     0.01     |
|   create_dof_constraints()         2         0.0898      0.044914    0.4939      0.246931    3.26     17.91    |
|   distribute_dofs()                2         0.0540      0.027023    0.2044      0.102203    1.96     7.41     |
|   dof_indices()                    7460      0.6932      0.000093    0.6932      0.000093    25.14    25.14    |
|   prepare_send_list()              2         0.0011      0.000535    0.0011      0.000535    0.04     0.04     |
|   reinit()                         2         0.1211      0.060551    0.1211      0.060551    4.39     4.39     |
|                                                                                                                |
| EquationSystems                                                                                                |
|   build_solution_vector()          1         0.0137      0.013744    0.0781      0.078150    0.50     2.83     |
|                                                                                                                |
| ExodusII_IO                                                                                                    |
|   write_nodal_data()               1         0.0190      0.018980    0.0190      0.018980    0.69     0.69     |
|                                                                                                                |
| FE                                                                                                             |
|   compute_shape_functions()        298       0.0066      0.000022    0.0066      0.000022    0.24     0.24     |
|   init_shape_functions()           84        0.0006      0.000008    0.0006      0.000008    0.02     0.02     |
|                                                                                                                |
| FEMap                                                                                                          |
|   compute_affine_map()             298       0.0056      0.000019    0.0056      0.000019    0.20     0.20     |
|   compute_face_map()               82        0.0116      0.000141    0.0116      0.000141    0.42     0.42     |
|   init_face_shape_functions()      1         0.0000      0.000049    0.0000      0.000049    0.00     0.00     |
|   init_reference_to_physical_map() 84        0.0132      0.000157    0.0132      0.000157    0.48     0.48     |
|                                                                                                                |
| Mesh                                                                                                           |
|   find_neighbors()                 1         0.0470      0.047031    0.0484      0.048425    1.71     1.76     |
|   renumber_nodes_and_elem()        2         0.0020      0.001012    0.0020      0.001012    0.07     0.07     |
|                                                                                                                |
| MeshCommunication                                                                                              |
|   compute_hilbert_indices()        2         0.0215      0.010772    0.0215      0.010772    0.78     0.78     |
|   find_global_indices()            2         0.0082      0.004123    0.0369      0.018431    0.30     1.34     |
|   parallel_sort()                  2         0.0036      0.001819    0.0056      0.002782    0.13     0.20     |
|                                                                                                                |
| MeshOutput                                                                                                     |
|   write_equation_systems()         1         0.0002      0.000157    0.0974      0.097398    0.01     3.53     |
|                                                                                                                |
| MeshTools::Generation                                                                                          |
|   build_cube()                     1         0.0103      0.010307    0.0103      0.010307    0.37     0.37     |
|                                                                                                                |
| MetisPartitioner                                                                                               |
|   partition()                      1         0.0989      0.098918    0.1164      0.116364    3.59     4.22     |
|                                                                                                                |
| Parallel                                                                                                       |
|   allgather()                      12        0.0185      0.001539    0.0188      0.001563    0.67     0.68     |
|   max(bool)                        1         0.0000      0.000008    0.0000      0.000008    0.00     0.00     |
|   max(scalar)                      137       0.0014      0.000010    0.0014      0.000010    0.05     0.05     |
|   max(vector)                      31        0.0004      0.000014    0.0011      0.000037    0.02     0.04     |
|   min(bool)                        157       0.0012      0.000008    0.0012      0.000008    0.04     0.04     |
|   min(scalar)                      128       0.1045      0.000816    0.1045      0.000816    3.79     3.79     |
|   min(vector)                      31        0.0005      0.000018    0.0016      0.000053    0.02     0.06     |
|   probe()                          176       0.0049      0.000028    0.0049      0.000028    0.18     0.18     |
|   receive()                        176       0.0012      0.000007    0.0062      0.000035    0.05     0.22     |
|   send()                           176       0.0006      0.000003    0.0006      0.000003    0.02     0.02     |
|   send_receive()                   180       0.0019      0.000010    0.0091      0.000050    0.07     0.33     |
|   sum()                            28        0.0041      0.000147    0.0062      0.000220    0.15     0.22     |
|                                                                                                                |
| Parallel::Request                                                                                              |
|   wait()                           176       0.0004      0.000002    0.0004      0.000002    0.01     0.01     |
|                                                                                                                |
| Partitioner                                                                                                    |
|   set_node_processor_ids()         1         0.0032      0.003160    0.0065      0.006536    0.11     0.24     |
|   set_parent_processor_ids()       1         0.0024      0.002420    0.0024      0.002420    0.09     0.09     |
|                                                                                                                |
| PetscLinearSolver                                                                                              |
|   solve()                          1         0.3740      0.374006    0.3740      0.374006    13.57    13.57    |
|                                                                                                                |
| System                                                                                                         |
|   assemble()                       1         0.9021      0.902138    1.0153      1.015295    32.72    36.83    |
 ----------------------------------------------------------------------------------------------------------------
| Totals:                            9961      2.7570                                          100.00            |
 ----------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example systems_of_equations_ex6:
*  mpirun -np 12 example-devel -ksp_type cg -pc_type bjacobi -sub_pc_type icc -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
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
