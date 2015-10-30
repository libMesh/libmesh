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
        #include "libmesh.h"
        #include "mesh.h"
        #include "mesh_generation.h"
        #include "exodusII_io.h"
        #include "gnuplot_io.h"
        #include "linear_implicit_system.h"
        #include "equation_systems.h"
        #include "fe.h"
        #include "quadrature_gauss.h"
        #include "dof_map.h"
        #include "sparse_matrix.h"
        #include "numeric_vector.h"
        #include "dense_matrix.h"
        #include "dense_submatrix.h"
        #include "dense_vector.h"
        #include "dense_subvector.h"
        #include "perf_log.h"
        #include "elem.h"
        #include "boundary_info.h"
        #include "zero_function.h"
        #include "dirichlet_boundaries.h"
        #include "string_to_enum.h"
        #include "getpot.h"
        
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
          libmesh_assert (system_name == "Elasticity");
          
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
        
          std::vector&lt;unsigned int&gt; dof_indices;
          std::vector&lt;unsigned int&gt; dof_indices_u;
          std::vector&lt;unsigned int&gt; dof_indices_v;
          std::vector&lt;unsigned int&gt; dof_indices_w;
        
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
                      boundary_id_type bc_id = mesh.boundary_info-&gt;boundary_id (elem,side);
                      if (bc_id==BoundaryInfo::invalid_id)
                          libmesh_error();
        
                      const std::vector&lt;std::vector&lt;Real&gt; &gt;&  phi_face = fe_face-&gt;get_phi();
                      const std::vector&lt;Real&gt;& JxW_face = fe_face-&gt;get_JxW();
        
                      fe_face-&gt;reinit(elem, side);
        
                      for (unsigned int qp=0; qp&lt;qface.n_points(); qp++)
                      {
                        if( bc_id == BOUNDARY_ID_MAX_X ) // Apply a traction
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
          std::vector&lt; std::vector&lt;unsigned int&gt; &gt; dof_indices_var(system.n_vars());
          std::vector&lt;unsigned int&gt; stress_dof_indices_var;
        
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
                unsigned int dof_index = stress_dof_indices_var[0];
                
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
            unsigned int dof_index = stress_dof_indices_var[0];
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
<br><br><br> <h1> The program without comments: </h1> 
<pre> 
  
  
  #include &lt;iostream&gt;
  #include &lt;algorithm&gt;
  #include &lt;math.h&gt;
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;mesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;mesh_generation.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;exodusII_io.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;gnuplot_io.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;linear_implicit_system.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;equation_systems.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;fe.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;quadrature_gauss.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;dof_map.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;sparse_matrix.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;numeric_vector.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;dense_matrix.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;dense_submatrix.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;dense_vector.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;dense_subvector.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;perf_log.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;elem.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;boundary_info.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;zero_function.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;dirichlet_boundaries.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;string_to_enum.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;getpot.h&quot;</FONT></B>
  
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
    libmesh_assert (system_name == <B><FONT COLOR="#BC8F8F">&quot;Elasticity&quot;</FONT></B>);
    
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
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B>&gt; dof_indices;
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B>&gt; dof_indices_u;
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B>&gt; dof_indices_v;
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B>&gt; dof_indices_w;
  
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
                boundary_id_type bc_id = mesh.boundary_info-&gt;boundary_id (elem,side);
                <B><FONT COLOR="#A020F0">if</FONT></B> (bc_id==BoundaryInfo::invalid_id)
                    libmesh_error();
  
                <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;Real&gt; &gt;&amp;  phi_face = fe_face-&gt;get_phi();
                <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Real&gt;&amp; JxW_face = fe_face-&gt;get_JxW();
  
                fe_face-&gt;reinit(elem, side);
  
                <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> qp=0; qp&lt;qface.n_points(); qp++)
                {
                  <B><FONT COLOR="#A020F0">if</FONT></B>( bc_id == BOUNDARY_ID_MAX_X ) <I><FONT COLOR="#B22222">// Apply a traction
</FONT></I>                  {
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
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt; std::vector&lt;<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B>&gt; &gt; dof_indices_var(system.n_vars());
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B>&gt; stress_dof_indices_var;
  
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
  
          <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> dof_index = stress_dof_indices_var[0];
          
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
      <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> dof_index = stress_dof_indices_var[0];
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
Linking systems_of_equations_ex6-opt...
***************************************************************
* Running  mpirun -np 6 ./systems_of_equations_ex6-opt
***************************************************************
 
 Mesh Information:
  mesh_dimension()=3
  spatial_dimension()=3
  n_nodes()=1845
    n_local_nodes()=348
  n_elem()=1280
    n_local_elem()=214
    n_active_elem()=1280
  n_subdomains()=1
  n_partitions()=6
  n_processors()=6
  n_threads()=1
  processor_id()=0

 EquationSystems
  n_systems()=2
   System #0, "Elasticity"
    Type "LinearImplicit"
    Variables="u" "v" "w" 
    Finite Element Types="LAGRANGE", "JACOBI_20_00" "LAGRANGE", "JACOBI_20_00" "LAGRANGE", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" "CARTESIAN" "CARTESIAN" 
    Approximation Orders="FIRST", "THIRD" "FIRST", "THIRD" "FIRST", "THIRD" 
    n_dofs()=5535
    n_local_dofs()=1044
    n_constrained_dofs()=135
    n_local_constrained_dofs()=135
    n_vectors()=1
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 59.5345
      Average Off-Processor Bandwidth <= 2.96552
      Maximum  On-Processor Bandwidth <= 81
      Maximum Off-Processor Bandwidth <= 36
    DofMap Constraints
      Number of DoF Constraints = 135
      Average DoF Constraint Length= 0
      Number of Node Constraints = 0
   System #1, "StressSystem"
    Type "Explicit"
    Variables="sigma_00" "sigma_01" "sigma_02" "sigma_10" "sigma_11" "sigma_12" "sigma_20" "sigma_21" "sigma_22" "vonMises" 
    Finite Element Types="MONOMIAL", "JACOBI_20_00" "MONOMIAL", "JACOBI_20_00" "MONOMIAL", "JACOBI_20_00" "MONOMIAL", "JACOBI_20_00" "MONOMIAL", "JACOBI_20_00" "MONOMIAL", "JACOBI_20_00" "MONOMIAL", "JACOBI_20_00" "MONOMIAL", "JACOBI_20_00" "MONOMIAL", "JACOBI_20_00" "MONOMIAL", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" "CARTESIAN" "CARTESIAN" "CARTESIAN" "CARTESIAN" "CARTESIAN" "CARTESIAN" "CARTESIAN" "CARTESIAN" "CARTESIAN" 
    Approximation Orders="CONSTANT", "THIRD" "CONSTANT", "THIRD" "CONSTANT", "THIRD" "CONSTANT", "THIRD" "CONSTANT", "THIRD" "CONSTANT", "THIRD" "CONSTANT", "THIRD" "CONSTANT", "THIRD" "CONSTANT", "THIRD" "CONSTANT", "THIRD" 
    n_dofs()=12800
    n_local_dofs()=2140
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


-------------------------------------------------------------------
| Processor id:   0                                                |
| Num Processors: 6                                                |
| Time:           Fri Aug 24 15:27:51 2012                         |
| OS:             Linux                                            |
| HostName:       daedalus                                         |
| OS Release:     2.6.32-34-generic                                |
| OS Version:     #76-Ubuntu SMP Tue Aug 30 17:05:01 UTC 2011      |
| Machine:        x86_64                                           |
| Username:       roystgnr                                         |
| Configuration:  ./configure run on Wed Aug 22 12:44:06 CDT 2012  |
-------------------------------------------------------------------
 ----------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=3.59305, Active time=3.40158                                                   |
 ----------------------------------------------------------------------------------------------------------------
| Event                              nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                              w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|----------------------------------------------------------------------------------------------------------------|
|                                                                                                                |
|                                                                                                                |
| DofMap                                                                                                         |
|   add_neighbors_to_send_list()     2         0.0004      0.000194    0.0005      0.000272    0.01     0.02     |
|   build_constraint_matrix()        214       0.0005      0.000002    0.0005      0.000002    0.01     0.01     |
|   build_sparsity()                 1         0.0078      0.007820    0.0088      0.008776    0.23     0.26     |
|   cnstrn_elem_mat_vec()            214       0.0077      0.000036    0.0077      0.000036    0.23     0.23     |
|   create_dof_constraints()         2         0.0055      0.002761    0.0070      0.003515    0.16     0.21     |
|   distribute_dofs()                2         0.0026      0.001313    0.0281      0.014069    0.08     0.83     |
|   dof_indices()                    11694     0.0045      0.000000    0.0045      0.000000    0.13     0.13     |
|   prepare_send_list()              2         0.0001      0.000029    0.0001      0.000029    0.00     0.00     |
|   reinit()                         2         0.0072      0.003589    0.0072      0.003589    0.21     0.21     |
|                                                                                                                |
| EquationSystems                                                                                                |
|   build_solution_vector()          1         0.0025      0.002518    0.1628      0.162818    0.07     4.79     |
|                                                                                                                |
| ExodusII_IO                                                                                                    |
|   write_nodal_data()               1         0.0028      0.002791    0.0028      0.002791    0.08     0.08     |
|                                                                                                                |
| FE                                                                                                             |
|   compute_shape_functions()        617       0.0010      0.000002    0.0010      0.000002    0.03     0.03     |
|   init_shape_functions()           191       0.0001      0.000001    0.0001      0.000001    0.00     0.00     |
|                                                                                                                |
| FEMap                                                                                                          |
|   compute_affine_map()             617       0.0008      0.000001    0.0008      0.000001    0.02     0.02     |
|   compute_face_map()               189       0.0003      0.000001    0.0003      0.000001    0.01     0.01     |
|   init_face_shape_functions()      1         0.0000      0.000010    0.0000      0.000010    0.00     0.00     |
|   init_reference_to_physical_map() 191       0.0010      0.000005    0.0010      0.000005    0.03     0.03     |
|                                                                                                                |
| Mesh                                                                                                           |
|   find_neighbors()                 1         0.0020      0.002036    0.0025      0.002550    0.06     0.07     |
|   renumber_nodes_and_elem()        2         0.0002      0.000082    0.0002      0.000082    0.00     0.00     |
|                                                                                                                |
| MeshCommunication                                                                                              |
|   compute_hilbert_indices()        2         0.0043      0.002129    0.0043      0.002129    0.13     0.13     |
|   find_global_indices()            2         0.0004      0.000224    0.0122      0.006120    0.01     0.36     |
|   parallel_sort()                  2         0.0011      0.000563    0.0065      0.003268    0.03     0.19     |
|                                                                                                                |
| MeshOutput                                                                                                     |
|   write_equation_systems()         1         0.0000      0.000034    0.1656      0.165643    0.00     4.87     |
|                                                                                                                |
| MeshTools::Generation                                                                                          |
|   build_cube()                     1         0.0007      0.000685    0.0007      0.000685    0.02     0.02     |
|                                                                                                                |
| MetisPartitioner                                                                                               |
|   partition()                      1         0.0028      0.002761    0.0082      0.008198    0.08     0.24     |
|                                                                                                                |
| Parallel                                                                                                       |
|   allgather()                      11        0.0149      0.001356    0.0149      0.001356    0.44     0.44     |
|   broadcast()                      1         0.0000      0.000012    0.0000      0.000012    0.00     0.00     |
|   gather()                         1         0.0000      0.000044    0.0000      0.000044    0.00     0.00     |
|   max(scalar)                      3         0.0080      0.002664    0.0080      0.002664    0.23     0.23     |
|   max(vector)                      2         0.0001      0.000051    0.0001      0.000051    0.00     0.00     |
|   min(vector)                      2         0.0006      0.000286    0.0006      0.000286    0.02     0.02     |
|   probe()                          70        0.0072      0.000103    0.0072      0.000103    0.21     0.21     |
|   receive()                        70        0.0002      0.000003    0.0074      0.000106    0.01     0.22     |
|   send()                           70        0.0002      0.000002    0.0002      0.000002    0.00     0.00     |
|   send_receive()                   74        0.0001      0.000002    0.0078      0.000105    0.00     0.23     |
|   sum()                            12        0.1651      0.013760    0.1651      0.013760    4.85     4.85     |
|                                                                                                                |
| Parallel::Request                                                                                              |
|   wait()                           70        0.0000      0.000001    0.0000      0.000001    0.00     0.00     |
|                                                                                                                |
| Partitioner                                                                                                    |
|   set_node_processor_ids()         1         0.0002      0.000199    0.0032      0.003240    0.01     0.10     |
|   set_parent_processor_ids()       1         0.0001      0.000077    0.0001      0.000077    0.00     0.00     |
|                                                                                                                |
| PetscLinearSolver                                                                                              |
|   solve()                          1         2.9961      2.996079    2.9961      2.996079    88.08    88.08    |
|                                                                                                                |
| System                                                                                                         |
|   assemble()                       1         0.1525      0.152523    0.1639      0.163899    4.48     4.82     |
 ----------------------------------------------------------------------------------------------------------------
| Totals:                            14343     3.4016                                          100.00            |
 ----------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running  mpirun -np 6 ./systems_of_equations_ex6-opt
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
