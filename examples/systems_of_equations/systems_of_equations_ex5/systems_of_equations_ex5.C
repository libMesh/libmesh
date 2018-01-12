// The libMesh Finite Element Library.
// Copyright (C) 2002-2018 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA



// <h1> Systems Example 5 - Linear Elastic Cantilever with Constraint </h1>
// \author David Knezevic
// \date 2012
//
// In this example we extend systems_of_equations_ex4 to enforce a constraint.
// We apply a uniform load on the top surface of the cantilever, and we
// determine the traction on the right boundary in order to obtain zero
// average vertical displacement on the right boundary of the domain.
//
// This constraint is enforced via a Lagrange multiplier (SCALAR variable).
// The system we solve, therefore, is of the form:
//   a(u,v) + \lambda g(v) = f(v)
//                    g(u) = 0
// Here \lambda tells us the traction required to satisfy the constraint.

// C++ include files that we need
#include <iostream>
#include <algorithm>
#include <math.h>

// libMesh includes
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

// Bring in everything from the libMesh namespace
using namespace libMesh;

// Matrix and right-hand side assemble
void assemble_elasticity(EquationSystems & es,
                         const std::string & system_name);

// Define the elasticity tensor, which is a fourth-order tensor
// i.e. it has four indices i, j, k, l
Real eval_elasticity_tensor(unsigned int i,
                            unsigned int j,
                            unsigned int k,
                            unsigned int l);

// Begin the main program.
int main (int argc, char ** argv)
{
  // Initialize libMesh and any dependent libraries
  LibMeshInit init (argc, argv);

  // This example requires a linear solver package.
  libmesh_example_requires(libMesh::default_solver_package() != INVALID_SOLVER_PACKAGE,
                           "--enable-petsc, --enable-trilinos, or --enable-eigen");

  // This example NaNs with the Eigen sparse linear solvers
  libmesh_example_requires(libMesh::default_solver_package() != EIGEN_SOLVERS, "--enable-petsc or --enable-laspack");

  // Initialize the cantilever mesh
  const unsigned int dim = 2;

  // Skip this 2D example if libMesh was compiled as 1D-only.
  libmesh_example_requires(dim <= LIBMESH_DIM, "2D support");

  // Create a 2D mesh distributed across the default MPI communicator.
  Mesh mesh(init.comm(), dim);
  MeshTools::Generation::build_square (mesh,
                                       50, 10,
                                       0., 1.,
                                       0., 0.2,
                                       QUAD9);

  // Print information about the mesh to the screen.
  mesh.print_info();

  // Create an equation systems object.
  EquationSystems equation_systems (mesh);

  // Declare the system and its variables.
  // Create a system named "Elasticity"
  LinearImplicitSystem & system =
    equation_systems.add_system<LinearImplicitSystem> ("Elasticity");

  // Add two displacement variables, u and v, to the system
  unsigned int u_var = system.add_variable("u", SECOND, LAGRANGE);
  unsigned int v_var = system.add_variable("v", SECOND, LAGRANGE);

  // Add a SCALAR variable for the Lagrange multiplier to enforce our constraint
  system.add_variable("lambda", FIRST, SCALAR);

  system.attach_assemble_function (assemble_elasticity);

  // Construct a Dirichlet boundary condition object
  // We impose a "clamped" boundary condition on the
  // "left" boundary, i.e. bc_id = 3
  std::set<boundary_id_type> boundary_ids;
  boundary_ids.insert(3);

  // Create a vector storing the variable numbers which the BC applies to
  std::vector<unsigned int> variables(2);
  variables[0] = u_var;
  variables[1] = v_var;

  // Create a ZeroFunction to initialize dirichlet_bc
  ZeroFunction<> zf;

  // Most DirichletBoundary users will want to supply a "locally
  // indexed" functor
  DirichletBoundary dirichlet_bc(boundary_ids, variables, zf,
                                 LOCAL_VARIABLE_ORDER);

  // We must add the Dirichlet boundary condition _before_
  // we call equation_systems.init()
  system.get_dof_map().add_dirichlet_boundary(dirichlet_bc);

  // Initialize the data structures for the equation system.
  equation_systems.init();

  // Print information about the system to the screen.
  equation_systems.print_info();

  // Solve the system
  system.solve();

  // Plot the solution
#ifdef LIBMESH_HAVE_EXODUS_API
  ExodusII_IO (mesh).write_equation_systems("displacement.e", equation_systems);
#endif // #ifdef LIBMESH_HAVE_EXODUS_API

  // All done.
  return 0;
}


void assemble_elasticity(EquationSystems & es,
                         const std::string & libmesh_dbg_var(system_name))
{
  libmesh_assert_equal_to (system_name, "Elasticity");

  const MeshBase & mesh = es.get_mesh();

  const unsigned int dim = mesh.mesh_dimension();

  LinearImplicitSystem & system = es.get_system<LinearImplicitSystem>("Elasticity");

  const unsigned int u_var = system.variable_number ("u");
  const unsigned int v_var = system.variable_number ("v");
  const unsigned int lambda_var = system.variable_number ("lambda");

  const DofMap & dof_map = system.get_dof_map();
  FEType fe_type = dof_map.variable_type(0);
  std::unique_ptr<FEBase> fe (FEBase::build(dim, fe_type));
  QGauss qrule (dim, fe_type.default_quadrature_order());
  fe->attach_quadrature_rule (&qrule);

  std::unique_ptr<FEBase> fe_face (FEBase::build(dim, fe_type));
  QGauss qface(dim-1, fe_type.default_quadrature_order());
  fe_face->attach_quadrature_rule (&qface);

  const std::vector<Real> & JxW = fe->get_JxW();
  const std::vector<std::vector<RealGradient>> & dphi = fe->get_dphi();

  DenseMatrix<Number> Ke;
  DenseVector<Number> Fe;

  DenseSubMatrix<Number>
    Kuu(Ke), Kuv(Ke),
    Kvu(Ke), Kvv(Ke);
  DenseSubMatrix<Number> Klambda_v(Ke), Kv_lambda(Ke);

  DenseSubVector<Number>
    Fu(Fe),
    Fv(Fe);

  std::vector<dof_id_type> dof_indices;
  std::vector<dof_id_type> dof_indices_u;
  std::vector<dof_id_type> dof_indices_v;
  std::vector<dof_id_type> dof_indices_lambda;

  for (const auto & elem : mesh.active_local_element_ptr_range())
    {
      dof_map.dof_indices (elem, dof_indices);
      dof_map.dof_indices (elem, dof_indices_u, u_var);
      dof_map.dof_indices (elem, dof_indices_v, v_var);
      dof_map.dof_indices (elem, dof_indices_lambda, lambda_var);

      const unsigned int n_dofs   = dof_indices.size();
      const unsigned int n_u_dofs = dof_indices_u.size();
      const unsigned int n_v_dofs = dof_indices_v.size();
      const unsigned int n_lambda_dofs = dof_indices_lambda.size();

      fe->reinit (elem);

      Ke.resize (n_dofs, n_dofs);
      Fe.resize (n_dofs);

      Kuu.reposition (u_var*n_u_dofs, u_var*n_u_dofs, n_u_dofs, n_u_dofs);
      Kuv.reposition (u_var*n_u_dofs, v_var*n_u_dofs, n_u_dofs, n_v_dofs);

      Kvu.reposition (v_var*n_v_dofs, u_var*n_v_dofs, n_v_dofs, n_u_dofs);
      Kvv.reposition (v_var*n_v_dofs, v_var*n_v_dofs, n_v_dofs, n_v_dofs);

      // Also, add a row and a column to enforce the constraint
      Kv_lambda.reposition (v_var*n_u_dofs, v_var*n_u_dofs+n_v_dofs, n_v_dofs, 1);
      Klambda_v.reposition (v_var*n_v_dofs+n_v_dofs, v_var*n_v_dofs, 1, n_v_dofs);

      Fu.reposition (u_var*n_u_dofs, n_u_dofs);
      Fv.reposition (v_var*n_u_dofs, n_v_dofs);

      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
        {
          for (unsigned int i=0; i<n_u_dofs; i++)
            for (unsigned int j=0; j<n_u_dofs; j++)
              {
                // Tensor indices
                unsigned int C_i, C_j, C_k, C_l;
                C_i=0, C_k=0;

                C_j=0, C_l=0;
                Kuu(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i, C_j, C_k, C_l) * dphi[i][qp](C_j)*dphi[j][qp](C_l));

                C_j=1, C_l=0;
                Kuu(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i, C_j, C_k, C_l) * dphi[i][qp](C_j)*dphi[j][qp](C_l));

                C_j=0, C_l=1;
                Kuu(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i, C_j, C_k, C_l) * dphi[i][qp](C_j)*dphi[j][qp](C_l));

                C_j=1, C_l=1;
                Kuu(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i, C_j, C_k, C_l) * dphi[i][qp](C_j)*dphi[j][qp](C_l));
              }

          for (unsigned int i=0; i<n_u_dofs; i++)
            for (unsigned int j=0; j<n_v_dofs; j++)
              {
                // Tensor indices
                unsigned int C_i, C_j, C_k, C_l;
                C_i=0, C_k=1;

                C_j=0, C_l=0;
                Kuv(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i, C_j, C_k, C_l) * dphi[i][qp](C_j)*dphi[j][qp](C_l));

                C_j=1, C_l=0;
                Kuv(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i, C_j, C_k, C_l) * dphi[i][qp](C_j)*dphi[j][qp](C_l));

                C_j=0, C_l=1;
                Kuv(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i, C_j, C_k, C_l) * dphi[i][qp](C_j)*dphi[j][qp](C_l));

                C_j=1, C_l=1;
                Kuv(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i, C_j, C_k, C_l) * dphi[i][qp](C_j)*dphi[j][qp](C_l));
              }

          for (unsigned int i=0; i<n_v_dofs; i++)
            for (unsigned int j=0; j<n_u_dofs; j++)
              {
                // Tensor indices
                unsigned int C_i, C_j, C_k, C_l;
                C_i=1, C_k=0;

                C_j=0, C_l=0;
                Kvu(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i, C_j, C_k, C_l) * dphi[i][qp](C_j)*dphi[j][qp](C_l));

                C_j=1, C_l=0;
                Kvu(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i, C_j, C_k, C_l) * dphi[i][qp](C_j)*dphi[j][qp](C_l));

                C_j=0, C_l=1;
                Kvu(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i, C_j, C_k, C_l) * dphi[i][qp](C_j)*dphi[j][qp](C_l));

                C_j=1, C_l=1;
                Kvu(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i, C_j, C_k, C_l) * dphi[i][qp](C_j)*dphi[j][qp](C_l));
              }

          for (unsigned int i=0; i<n_v_dofs; i++)
            for (unsigned int j=0; j<n_v_dofs; j++)
              {
                // Tensor indices
                unsigned int C_i, C_j, C_k, C_l;
                C_i=1, C_k=1;

                C_j=0, C_l=0;
                Kvv(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i, C_j, C_k, C_l) * dphi[i][qp](C_j)*dphi[j][qp](C_l));

                C_j=1, C_l=0;
                Kvv(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i, C_j, C_k, C_l) * dphi[i][qp](C_j)*dphi[j][qp](C_l));

                C_j=0, C_l=1;
                Kvv(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i, C_j, C_k, C_l) * dphi[i][qp](C_j)*dphi[j][qp](C_l));

                C_j=1, C_l=1;
                Kvv(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i, C_j, C_k, C_l) * dphi[i][qp](C_j)*dphi[j][qp](C_l));
              }
        }

      {
        std::vector<boundary_id_type> bc_ids;
        for (auto side : elem->side_index_range())
          if (elem->neighbor_ptr(side) == libmesh_nullptr)
            {
              mesh.get_boundary_info().boundary_ids (elem, side, bc_ids);

              const std::vector<std::vector<Real>> & phi_face = fe_face->get_phi();
              const std::vector<Real> & JxW_face = fe_face->get_JxW();

              fe_face->reinit(elem, side);

              for (std::vector<boundary_id_type>::const_iterator b =
                     bc_ids.begin(); b != bc_ids.end(); ++b)
                {
                  const boundary_id_type bc_id = *b;
                  for (unsigned int qp=0; qp<qface.n_points(); qp++)
                    {
                      // Add the loading
                      if (bc_id == 2)
                        for (unsigned int i=0; i<n_v_dofs; i++)
                          Fv(i) += JxW_face[qp] * (-1.) * phi_face[i][qp];

                      // Add the constraint contributions
                      if (bc_id == 1)
                        {
                          for (unsigned int i=0; i<n_v_dofs; i++)
                            for (unsigned int j=0; j<n_lambda_dofs; j++)
                              Kv_lambda(i,j) += JxW_face[qp] * (-1.) * phi_face[i][qp];

                          for (unsigned int i=0; i<n_lambda_dofs; i++)
                            for (unsigned int j=0; j<n_v_dofs; j++)
                              Klambda_v(i,j) += JxW_face[qp] * (-1.) * phi_face[j][qp];
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

Real eval_elasticity_tensor(unsigned int i,
                            unsigned int j,
                            unsigned int k,
                            unsigned int l)
{
  // Define the Poisson ratio
  const Real nu = 0.3;

  // Define the Lame constants (lambda_1 and lambda_2) based on Poisson ratio
  const Real lambda_1 = nu / ((1. + nu) * (1. - 2.*nu));
  const Real lambda_2 = 0.5 / (1 + nu);

  // Define the Kronecker delta functions that we need here
  Real delta_ij = (i == j) ? 1. : 0.;
  Real delta_il = (i == l) ? 1. : 0.;
  Real delta_ik = (i == k) ? 1. : 0.;
  Real delta_jl = (j == l) ? 1. : 0.;
  Real delta_jk = (j == k) ? 1. : 0.;
  Real delta_kl = (k == l) ? 1. : 0.;

  return lambda_1 * delta_ij * delta_kl + lambda_2 * (delta_ik * delta_jl + delta_il * delta_jk);
}
