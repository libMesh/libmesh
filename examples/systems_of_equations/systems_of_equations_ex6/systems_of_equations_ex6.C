/* The libMesh Finite Element Library. */
/* Copyright (C) 2003  Benjamin S. Kirk */

/* This library is free software; you can redistribute it and/or */
/* modify it under the terms of the GNU Lesser General Public */
/* License as published by the Free Software Foundation; either */
/* version 2.1 of the License, or (at your option) any later version. */

/* This library is distributed in the hope that it will be useful, */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU */
/* Lesser General Public License for more details. */

/* You should have received a copy of the GNU Lesser General Public */
/* License along with this library; if not, write to the Free Software */
/* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA */



 // <h1> Systems Example 6 - 3D Linear Elastic Cantilever </h1>
 //      By David Knezevic
 //
 // This is a 3D version of systems_of_equations_ex4.
 // In this case we also compute and plot stresses.


// C++ include files that we need
#include <iostream>
#include <algorithm>
#include <math.h>

// libMesh includes
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

// boundary IDs
#define BOUNDARY_ID_MIN_Z 0
#define BOUNDARY_ID_MIN_Y 1
#define BOUNDARY_ID_MAX_X 2
#define BOUNDARY_ID_MAX_Y 3
#define BOUNDARY_ID_MIN_X 4
#define BOUNDARY_ID_MAX_Z 5
#define NODE_BOUNDARY_ID 10
#define EDGE_BOUNDARY_ID 20

// Bring in everything from the libMesh namespace
using namespace libMesh;

// Matrix and right-hand side assembly
void assemble_elasticity(EquationSystems& es,
                         const std::string& system_name);

// Post-process the solution to compute stresses
void compute_stresses(EquationSystems& es);

// The Kronecker delta function, used in eval_elasticity_tensor
Real kronecker_delta(unsigned int i,
                     unsigned int j);

// Define the elasticity tensor, which is a fourth-order tensor
// i.e. it has four indices i,j,k,l
Real eval_elasticity_tensor(unsigned int i,
                            unsigned int j,
                            unsigned int k,
                            unsigned int l);

// Begin the main program.
int main (int argc, char** argv)
{
  // Initialize libMesh and any dependent libaries
  LibMeshInit init (argc, argv);

  // Initialize the cantilever mesh
  const unsigned int dim = 3;

  // Make sure libMesh was compiled for 3D
  libmesh_example_assert(dim == LIBMESH_DIM, "3D support");

  // Create a 3D mesh distributed across the default MPI communicator.
  Mesh mesh(init.comm(), dim);
  MeshTools::Generation::build_cube (mesh,
                                     40,
                                     10,
                                     5,
                                     0., 1.*x_scaling,
                                     0., 0.3,
                                     0., 0.1,
                                     HEX8);


  // Print information about the mesh to the screen.
  mesh.print_info();

  // Let's add some node and edge boundary conditions
  MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
  for ( ; el != end_el; ++el)
  {
    const Elem* elem = *el;
    
    unsigned int side_max_x = 0, side_min_y = 0,
                 side_max_y = 0, side_max_z = 0;
    bool found_side_max_x = false, found_side_max_y = false,
         found_side_min_y = false, found_side_max_z = false;
    for(unsigned int side=0; side<elem->n_sides(); side++)
    {
      if( mesh.boundary_info->has_boundary_id(elem, side, BOUNDARY_ID_MAX_X))
      {
        side_max_x = side;
        found_side_max_x = true;
      }

      if( mesh.boundary_info->has_boundary_id(elem, side, BOUNDARY_ID_MIN_Y))
      {
        side_min_y = side;
        found_side_min_y = true;
      }

      if( mesh.boundary_info->has_boundary_id(elem, side, BOUNDARY_ID_MAX_Y))
      {
        side_max_y = side;
        found_side_max_y = true;
      }
      
      if( mesh.boundary_info->has_boundary_id(elem, side, BOUNDARY_ID_MAX_Z))
      {
        side_max_z = side;
        found_side_max_z = true;
      }
    }
    
    // If elem has sides on boundaries
    // BOUNDARY_ID_MAX_X, BOUNDARY_ID_MAX_Y, BOUNDARY_ID_MAX_Z
    // then let's set a node boundary condition
    if(found_side_max_x && found_side_max_y && found_side_max_z)
    {
      for(unsigned int n=0; n<elem->n_nodes(); n++)
      {
        if (elem->is_node_on_side(n, side_max_x) &&
            elem->is_node_on_side(n, side_max_y) &&
            elem->is_node_on_side(n, side_max_z) )
        {
          mesh.boundary_info->add_node(elem->get_node(n), NODE_BOUNDARY_ID);
        }
      }
    }
    
    
    // If elem has sides on boundaries
    // BOUNDARY_ID_MAX_X and BOUNDARY_ID_MIN_Y
    // then let's set an edge boundary condition
    if(found_side_max_x && found_side_min_y)
    {
      for(unsigned int e=0; e<elem->n_edges(); e++)
      {
        if (elem->is_edge_on_side(e, side_max_x) &&
            elem->is_edge_on_side(e, side_min_y) )
        {
          mesh.boundary_info->add_edge(elem, e, EDGE_BOUNDARY_ID);
        }
      }
    }
  }

  // Create an equation systems object.
  EquationSystems equation_systems (mesh);

  // Declare the system and its variables.
  // Create a system named "Elasticity"
  LinearImplicitSystem& system =
    equation_systems.add_system<LinearImplicitSystem> ("Elasticity");

  // Add three displacement variables, u and v, to the system
  unsigned int u_var = system.add_variable("u", FIRST, LAGRANGE);
  unsigned int v_var = system.add_variable("v", FIRST, LAGRANGE);
  unsigned int w_var = system.add_variable("w", FIRST, LAGRANGE);

  system.attach_assemble_function (assemble_elasticity);


  std::set<boundary_id_type> boundary_ids;
  boundary_ids.insert(BOUNDARY_ID_MIN_X);
  boundary_ids.insert(NODE_BOUNDARY_ID);
  boundary_ids.insert(EDGE_BOUNDARY_ID);

  // Create a vector storing the variable numbers which the BC applies to
  std::vector<unsigned int> variables;
  variables.push_back(u_var);
  variables.push_back(v_var);
  variables.push_back(w_var);

  // Create a ZeroFunction to initialize dirichlet_bc
  ZeroFunction<> zf;

  DirichletBoundary dirichlet_bc(boundary_ids,
                                 variables,
                                 &zf);

  // We must add the Dirichlet boundary condition _before_
  // we call equation_systems.init()
  system.get_dof_map().add_dirichlet_boundary(dirichlet_bc);

  // Also, initialize an ExplicitSystem to store stresses
  ExplicitSystem& stress_system =
    equation_systems.add_system<ExplicitSystem> ("StressSystem");

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

  // Initialize the data structures for the equation system.
  equation_systems.init();

  // Print information about the system to the screen.
  equation_systems.print_info();

  // Solve the system
  system.solve();

  // Post-process the solution to compute the stresses
  compute_stresses(equation_systems);

  // Plot the solution
#ifdef LIBMESH_HAVE_EXODUS_API
  ExodusII_IO (mesh).write_equation_systems("displacement.e",equation_systems);
#endif // #ifdef LIBMESH_HAVE_EXODUS_API

  // All done.
  return 0;
}


void assemble_elasticity(EquationSystems& es,
                         const std::string& system_name)
{
  libmesh_assert_equal_to (system_name, "Elasticity");

  const MeshBase& mesh = es.get_mesh();

  const unsigned int dim = mesh.mesh_dimension();

  LinearImplicitSystem& system = es.get_system<LinearImplicitSystem>("Elasticity");

  const unsigned int n_components = 3;
  const unsigned int u_var = system.variable_number ("u");
  const unsigned int v_var = system.variable_number ("v");
  const unsigned int w_var = system.variable_number ("w");

  const DofMap& dof_map = system.get_dof_map();
  FEType fe_type = dof_map.variable_type(u_var);
  AutoPtr<FEBase> fe (FEBase::build(dim, fe_type));
  QGauss qrule (dim, fe_type.default_quadrature_order());
  fe->attach_quadrature_rule (&qrule);

  AutoPtr<FEBase> fe_face (FEBase::build(dim, fe_type));
  QGauss qface(dim-1, fe_type.default_quadrature_order());
  fe_face->attach_quadrature_rule (&qface);

  const std::vector<Real>& JxW = fe->get_JxW();
  const std::vector<std::vector<Real> >& phi = fe->get_phi();
  const std::vector<std::vector<RealGradient> >& dphi = fe->get_dphi();

  DenseMatrix<Number> Ke;
  DenseVector<Number> Fe;

  DenseSubMatrix<Number>
    Kuu(Ke), Kuv(Ke), Kuw(Ke),
    Kvu(Ke), Kvv(Ke), Kvw(Ke),
    Kwu(Ke), Kwv(Ke), Kww(Ke);

  DenseSubVector<Number>
    Fu(Fe),
    Fv(Fe),
    Fw(Fe);

  std::vector<dof_id_type> dof_indices;
  std::vector<dof_id_type> dof_indices_u;
  std::vector<dof_id_type> dof_indices_v;
  std::vector<dof_id_type> dof_indices_w;

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

      fe->reinit (elem);

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

      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
      {
          for (unsigned int i=0; i<n_u_dofs; i++)
            for (unsigned int j=0; j<n_u_dofs; j++)
            {
              // Tensor indices
              unsigned int C_i=0, C_k=0;

              for(unsigned int C_j=0; C_j<n_components; C_j++)
                for(unsigned int C_l=0; C_l<n_components; C_l++)
                {
                  Kuu(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i,C_j,C_k,C_l) * dphi[i][qp](C_j)*dphi[j][qp](C_l));
                }
            }

          for (unsigned int i=0; i<n_u_dofs; i++)
            for (unsigned int j=0; j<n_v_dofs; j++)
            {
              // Tensor indices
              unsigned int C_i=0, C_k=1;

              for(unsigned int C_j=0; C_j<n_components; C_j++)
                for(unsigned int C_l=0; C_l<n_components; C_l++)
                {
                  Kuv(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i,C_j,C_k,C_l) * dphi[i][qp](C_j)*dphi[j][qp](C_l));
                }
            }

          for (unsigned int i=0; i<n_u_dofs; i++)
            for (unsigned int j=0; j<n_w_dofs; j++)
            {
              // Tensor indices
              unsigned int C_i=0, C_k=2;

              for(unsigned int C_j=0; C_j<n_components; C_j++)
                for(unsigned int C_l=0; C_l<n_components; C_l++)
                {
                  Kuw(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i,C_j,C_k,C_l) * dphi[i][qp](C_j)*dphi[j][qp](C_l));
                }
            }

          for (unsigned int i=0; i<n_v_dofs; i++)
            for (unsigned int j=0; j<n_u_dofs; j++)
            {
              // Tensor indices
              unsigned int C_i = 1, C_k = 0;

              for(unsigned int C_j=0; C_j<n_components; C_j++)
                for(unsigned int C_l=0; C_l<n_components; C_l++)
                {
                  Kvu(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i,C_j,C_k,C_l) * dphi[i][qp](C_j)*dphi[j][qp](C_l));
                }
            }

          for (unsigned int i=0; i<n_v_dofs; i++)
            for (unsigned int j=0; j<n_v_dofs; j++)
            {
              // Tensor indices
              unsigned int C_i = 1, C_k = 1;

              for(unsigned int C_j=0; C_j<n_components; C_j++)
                for(unsigned int C_l=0; C_l<n_components; C_l++)
                {
                  Kvv(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i,C_j,C_k,C_l) * dphi[i][qp](C_j)*dphi[j][qp](C_l));
                }
            }

          for (unsigned int i=0; i<n_v_dofs; i++)
            for (unsigned int j=0; j<n_w_dofs; j++)
            {
              // Tensor indices
              unsigned int C_i = 1, C_k = 2;

              for(unsigned int C_j=0; C_j<n_components; C_j++)
                for(unsigned int C_l=0; C_l<n_components; C_l++)
                {
                  Kvw(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i,C_j,C_k,C_l) * dphi[i][qp](C_j)*dphi[j][qp](C_l));
                }
            }

          for (unsigned int i=0; i<n_w_dofs; i++)
            for (unsigned int j=0; j<n_u_dofs; j++)
            {
              // Tensor indices
              unsigned int C_i = 2, C_k = 0;

              for(unsigned int C_j=0; C_j<n_components; C_j++)
                for(unsigned int C_l=0; C_l<n_components; C_l++)
                {
                  Kwu(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i,C_j,C_k,C_l) * dphi[i][qp](C_j)*dphi[j][qp](C_l));
                }
            }

          for (unsigned int i=0; i<n_w_dofs; i++)
            for (unsigned int j=0; j<n_v_dofs; j++)
            {
              // Tensor indices
              unsigned int C_i = 2, C_k = 1;

              for(unsigned int C_j=0; C_j<n_components; C_j++)
                for(unsigned int C_l=0; C_l<n_components; C_l++)
                {
                  Kwv(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i,C_j,C_k,C_l) * dphi[i][qp](C_j)*dphi[j][qp](C_l));
                }
            }

          for (unsigned int i=0; i<n_w_dofs; i++)
            for (unsigned int j=0; j<n_w_dofs; j++)
            {
              // Tensor indices
              unsigned int C_i = 2, C_k = 2;

              for(unsigned int C_j=0; C_j<n_components; C_j++)
                for(unsigned int C_l=0; C_l<n_components; C_l++)
                {
                  Kww(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i,C_j,C_k,C_l) * dphi[i][qp](C_j)*dphi[j][qp](C_l));
                }
            }
            
          // Volumetric load
          for (unsigned int i=0; i<n_w_dofs; i++)
          {
            Fw(i) -= JxW[qp] * phi[i][qp];
          }
      }

      {
        for (unsigned int side=0; side<elem->n_sides(); side++)
          if (elem->neighbor(side) == NULL)
            {
              const std::vector<std::vector<Real> >&  phi_face = fe_face->get_phi();
              const std::vector<Real>& JxW_face = fe_face->get_JxW();

              fe_face->reinit(elem, side);

              for (unsigned int qp=0; qp<qface.n_points(); qp++)
              {
                // Apply a traction
                if( mesh.boundary_info->has_boundary_id(elem, side, BOUNDARY_ID_MAX_X) )
                {
                  for (unsigned int i=0; i<n_v_dofs; i++)
                  {
                    Fu(i) += JxW_face[qp] * x_load * phi_face[i][qp];
                  }
                  for (unsigned int i=0; i<n_v_dofs; i++)
                  {
                    Fv(i) += JxW_face[qp] * y_load * phi_face[i][qp];
                  }
                  for (unsigned int i=0; i<n_v_dofs; i++)
                  {
                    Fw(i) += JxW_face[qp] * z_load * phi_face[i][qp];
                  }
                }
              }
            }
      }

      dof_map.constrain_element_matrix_and_vector (Ke, Fe, dof_indices);

      system.matrix->add_matrix (Ke, dof_indices);
      system.rhs->add_vector    (Fe, dof_indices);
    }
}

void compute_stresses(EquationSystems& es)
{
  const MeshBase& mesh = es.get_mesh();

  const unsigned int dim = mesh.mesh_dimension();

  LinearImplicitSystem& system = es.get_system<LinearImplicitSystem>("Elasticity");

  unsigned int displacement_vars[3];
  displacement_vars[0] = system.variable_number ("u");
  displacement_vars[1] = system.variable_number ("v");
  displacement_vars[2] = system.variable_number ("w");
  const unsigned int u_var = system.variable_number ("u");

  const DofMap& dof_map = system.get_dof_map();
  FEType fe_type = dof_map.variable_type(u_var);
  AutoPtr<FEBase> fe (FEBase::build(dim, fe_type));
  QGauss qrule (dim, fe_type.default_quadrature_order());
  fe->attach_quadrature_rule (&qrule);

  const std::vector<Real>& JxW = fe->get_JxW();
  const std::vector<std::vector<RealGradient> >& dphi = fe->get_dphi();

  // Also, get a reference to the ExplicitSystem
  ExplicitSystem& stress_system = es.get_system<ExplicitSystem>("StressSystem");
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

  // Storage for the stress dof indices on each element
  std::vector< std::vector<dof_id_type> > dof_indices_var(system.n_vars());
  std::vector<dof_id_type> stress_dof_indices_var;

  // To store the stress tensor on each element
  DenseMatrix<Number> elem_sigma;

  MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

  for ( ; el != end_el; ++el)
  {
    const Elem* elem = *el;

    for(unsigned int var=0; var<3; var++)
    {
      dof_map.dof_indices (elem, dof_indices_var[var], displacement_vars[var]);
    }

    fe->reinit (elem);

    elem_sigma.resize(3,3);

    for (unsigned int qp=0; qp<qrule.n_points(); qp++)
    {
      for(unsigned int C_i=0; C_i<3; C_i++)
        for(unsigned int C_j=0; C_j<3; C_j++)
          for(unsigned int C_k=0; C_k<3; C_k++)
          {
            const unsigned int n_var_dofs = dof_indices_var[C_k].size();

            // Get the gradient at this quadrature point
            Gradient displacement_gradient;
            for(unsigned int l=0; l<n_var_dofs; l++)
            {
              displacement_gradient.add_scaled(dphi[l][qp], system.current_solution(dof_indices_var[C_k][l]));
            }

            for(unsigned int C_l=0; C_l<3; C_l++)
            {
              elem_sigma(C_i,C_j) += JxW[qp]*( eval_elasticity_tensor(C_i,C_j,C_k,C_l) * displacement_gradient(C_l) );
            }

          }
    }

    // Get the average stresses by dividing by the element volume
    elem_sigma.scale(1./elem->volume());

    // load elem_sigma data into stress_system
    for(unsigned int i=0; i<3; i++)
      for(unsigned int j=0; j<3; j++)
      {
        stress_dof_map.dof_indices (elem, stress_dof_indices_var, sigma_vars[i][j]);

        // We are using CONSTANT MONOMIAL basis functions, hence we only need to get
        // one dof index per variable
        dof_id_type dof_index = stress_dof_indices_var[0];

        if( (stress_system.solution->first_local_index() <= dof_index) &&
            (dof_index < stress_system.solution->last_local_index()) )
        {
          stress_system.solution->set(dof_index, elem_sigma(i,j));
        }

      }

    // Also, the von Mises stress
    Number vonMises_value = std::sqrt( 0.5*( pow(elem_sigma(0,0) - elem_sigma(1,1),2.) +
                                             pow(elem_sigma(1,1) - elem_sigma(2,2),2.) +
                                             pow(elem_sigma(2,2) - elem_sigma(0,0),2.) +
                                             6.*(pow(elem_sigma(0,1),2.) + pow(elem_sigma(1,2),2.) + pow(elem_sigma(2,0),2.))
                                           ) );
    stress_dof_map.dof_indices (elem, stress_dof_indices_var, vonMises_var);
    dof_id_type dof_index = stress_dof_indices_var[0];
    if( (stress_system.solution->first_local_index() <= dof_index) &&
        (dof_index < stress_system.solution->last_local_index()) )
    {
      stress_system.solution->set(dof_index, vonMises_value);
    }

  }

  // Should call close and update when we set vector entries directly
  stress_system.solution->close();
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
  // Define the Poisson ratio and Young's modulus
  const Real nu = 0.3;
  const Real E  = 1.;

  // Define the Lame constants (lambda_1 and lambda_2) based on nu and E
  const Real lambda_1 = E * nu / ( (1. + nu) * (1. - 2.*nu) );
  const Real lambda_2 = 0.5 * E / (1. + nu);

  return lambda_1 * kronecker_delta(i,j) * kronecker_delta(k,l)
       + lambda_2 * (kronecker_delta(i,k) * kronecker_delta(j,l) + kronecker_delta(i,l) * kronecker_delta(j,k));
}
