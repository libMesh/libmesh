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



// <h1>Systems Example 1 - Stokes Equations</h1>
// \author Benjamin S. Kirk
// \date 2003
//
// This example shows how a simple, linear system of equations
// can be solved in parallel.  The system of equations are the familiar
// Stokes equations for low-speed incompressible fluid flow.

// C++ include files that we need
#include <iostream>
#include <algorithm>
#include <math.h>

// Basic include file needed for the mesh functionality.
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/equation_systems.h"
#include "libmesh/fe.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/dof_map.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/linear_implicit_system.h"

// For systems of equations the DenseSubMatrix
// and DenseSubVector provide convenient ways for
// assembling the element matrix and vector on a
// component-by-component basis.
#include "libmesh/dense_submatrix.h"
#include "libmesh/dense_subvector.h"

// The definition of a geometric element
#include "libmesh/elem.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;

// Function prototype.  This function will assemble the system
// matrix and right-hand-side.
void assemble_stokes (EquationSystems & es,
                      const std::string & system_name);

// The main program.
int main (int argc, char ** argv)
{
  // Initialize libMesh.
  LibMeshInit init (argc, argv);

  // This example requires a linear solver package.
  libmesh_example_requires(libMesh::default_solver_package() != INVALID_SOLVER_PACKAGE,
                           "--enable-petsc, --enable-trilinos, or --enable-eigen");

  // Skip this 2D example if libMesh was compiled as 1D-only.
  libmesh_example_requires(2 <= LIBMESH_DIM, "2D support");

  // This example NaNs with the Eigen sparse linear solvers and
  // Trilinos solvers, but should work OK with either PETSc or
  // Laspack.
  libmesh_example_requires(libMesh::default_solver_package() != EIGEN_SOLVERS, "--enable-petsc or --enable-laspack");
  libmesh_example_requires(libMesh::default_solver_package() != TRILINOS_SOLVERS, "--enable-petsc or --enable-laspack");

  // Create a mesh, with dimension to be overridden later, distributed
  // across the default MPI communicator.
  Mesh mesh(init.comm());

  // Use the MeshTools::Generation mesh generator to create a uniform
  // 2D grid on the square [-1,1]^2.  We instruct the mesh generator
  // to build a mesh of 8x8 Quad9 elements.  Building these
  // higher-order elements allows us to use higher-order
  // approximation, as in example 3.
  MeshTools::Generation::build_square (mesh,
                                       15, 15,
                                       0., 1.,
                                       0., 1.,
                                       QUAD9);

  // Print information about the mesh to the screen.
  mesh.print_info();

  // Create an equation systems object.
  EquationSystems equation_systems (mesh);

  // Declare the system and its variables.
  // Create a transient system named "Stokes"
  LinearImplicitSystem & system =
    equation_systems.add_system<LinearImplicitSystem> ("Stokes");

  // Add the variables "u" & "v" to "Stokes".  They
  // will be approximated using second-order approximation.
  system.add_variable ("u", SECOND);
  system.add_variable ("v", SECOND);

  // Add the variable "p" to "Stokes". This will
  // be approximated with a first-order basis,
  // providing an LBB-stable pressure-velocity pair.
  system.add_variable ("p", FIRST);

  // Give the system a pointer to the matrix assembly
  // function.
  system.attach_assemble_function (assemble_stokes);

  // Initialize the data structures for the equation system.
  equation_systems.init ();

  equation_systems.parameters.set<unsigned int>("linear solver maximum iterations") = 250;
  equation_systems.parameters.set<Real>        ("linear solver tolerance") = TOLERANCE;

  // Prints information about the system to the screen.
  equation_systems.print_info();

  // Assemble & solve the linear system,
  // then write the solution.
  equation_systems.get_system("Stokes").solve();

#ifdef LIBMESH_HAVE_EXODUS_API
  ExodusII_IO(mesh).write_equation_systems ("out.e",
                                            equation_systems);
#endif // #ifdef LIBMESH_HAVE_EXODUS_API

  // All done.
  return 0;
}

void assemble_stokes (EquationSystems & es,
                      const std::string & libmesh_dbg_var(system_name))
{
  // It is a good idea to make sure we are assembling
  // the proper system.
  libmesh_assert_equal_to (system_name, "Stokes");

  // Get a constant reference to the mesh object.
  const MeshBase & mesh = es.get_mesh();

  // The dimension that we are running
  const unsigned int dim = mesh.mesh_dimension();

  // Get a reference to the Convection-Diffusion system object.
  LinearImplicitSystem & system =
    es.get_system<LinearImplicitSystem> ("Stokes");

  // Numeric ids corresponding to each variable in the system
  const unsigned int u_var = system.variable_number ("u");
  const unsigned int v_var = system.variable_number ("v");
  const unsigned int p_var = system.variable_number ("p");

  // Get the Finite Element type for "u".  Note this will be
  // the same as the type for "v".
  FEType fe_vel_type = system.variable_type(u_var);

  // Get the Finite Element type for "p".
  FEType fe_pres_type = system.variable_type(p_var);

  // Build a Finite Element object of the specified type for
  // the velocity variables.
  std::unique_ptr<FEBase> fe_vel  (FEBase::build(dim, fe_vel_type));

  // Build a Finite Element object of the specified type for
  // the pressure variables.
  std::unique_ptr<FEBase> fe_pres (FEBase::build(dim, fe_pres_type));

  // A Gauss quadrature rule for numerical integration.
  // Let the FEType object decide what order rule is appropriate.
  QGauss qrule (dim, fe_vel_type.default_quadrature_order());

  // Tell the finite element objects to use our quadrature rule.
  fe_vel->attach_quadrature_rule (&qrule);
  fe_pres->attach_quadrature_rule (&qrule);

  // Here we define some references to cell-specific data that
  // will be used to assemble the linear system.
  //
  // The element Jacobian * quadrature weight at each integration point.
  const std::vector<Real> & JxW = fe_vel->get_JxW();

  // The element shape function gradients for the velocity
  // variables evaluated at the quadrature points.
  const std::vector<std::vector<RealGradient>> & dphi = fe_vel->get_dphi();

  // The element shape functions for the pressure variable
  // evaluated at the quadrature points.
  const std::vector<std::vector<Real>> & psi = fe_pres->get_phi();

  // A reference to the DofMap object for this system.  The DofMap
  // object handles the index translation from node and element numbers
  // to degree of freedom numbers.  We will talk more about the DofMap
  // in future examples.
  const DofMap & dof_map = system.get_dof_map();

  // Define data structures to contain the element matrix
  // and right-hand-side vector contribution.  Following
  // basic finite element terminology we will denote these
  // "Ke" and "Fe".
  DenseMatrix<Number> Ke;
  DenseVector<Number> Fe;

  DenseSubMatrix<Number>
    Kuu(Ke), Kuv(Ke), Kup(Ke),
    Kvu(Ke), Kvv(Ke), Kvp(Ke),
    Kpu(Ke), Kpv(Ke), Kpp(Ke);

  DenseSubVector<Number>
    Fu(Fe),
    Fv(Fe),
    Fp(Fe);

  // This vector will hold the degree of freedom indices for
  // the element.  These define where in the global system
  // the element degrees of freedom get mapped.
  std::vector<dof_id_type> dof_indices;
  std::vector<dof_id_type> dof_indices_u;
  std::vector<dof_id_type> dof_indices_v;
  std::vector<dof_id_type> dof_indices_p;

  // Now we will loop over all the elements in the mesh that
  // live on the local processor. We will compute the element
  // matrix and right-hand-side contribution.  In case users later
  // modify this program to include refinement, we will be safe and
  // will only consider the active elements; hence we use a variant of
  // the active_elem_iterator.
  for (const auto & elem : mesh.active_local_element_ptr_range())
    {
      // Get the degree of freedom indices for the
      // current element.  These define where in the global
      // matrix and right-hand-side this element will
      // contribute to.
      dof_map.dof_indices (elem, dof_indices);
      dof_map.dof_indices (elem, dof_indices_u, u_var);
      dof_map.dof_indices (elem, dof_indices_v, v_var);
      dof_map.dof_indices (elem, dof_indices_p, p_var);

      const unsigned int n_dofs   = dof_indices.size();
      const unsigned int n_u_dofs = dof_indices_u.size();
      const unsigned int n_v_dofs = dof_indices_v.size();
      const unsigned int n_p_dofs = dof_indices_p.size();

      // Compute the element-specific data for the current
      // element.  This involves computing the location of the
      // quadrature points (q_point) and the shape functions
      // (phi, dphi) for the current element.
      fe_vel->reinit  (elem);
      fe_pres->reinit (elem);

      // Zero the element matrix and right-hand side before
      // summing them.  We use the resize member here because
      // the number of degrees of freedom might have changed from
      // the last element.  Note that this will be the case if the
      // element type is different (i.e. the last element was a
      // triangle, now we are on a quadrilateral).
      Ke.resize (n_dofs, n_dofs);
      Fe.resize (n_dofs);

      // Reposition the submatrices...  The idea is this:
      //
      //         -           -          -  -
      //        | Kuu Kuv Kup |        | Fu |
      //   Ke = | Kvu Kvv Kvp |;  Fe = | Fv |
      //        | Kpu Kpv Kpp |        | Fp |
      //         -           -          -  -
      //
      // The DenseSubMatrix.reposition () member takes the
      // (row_offset, column_offset, row_size, column_size).
      //
      // Similarly, the DenseSubVector.reposition () member
      // takes the (row_offset, row_size)
      Kuu.reposition (u_var*n_u_dofs, u_var*n_u_dofs, n_u_dofs, n_u_dofs);
      Kuv.reposition (u_var*n_u_dofs, v_var*n_u_dofs, n_u_dofs, n_v_dofs);
      Kup.reposition (u_var*n_u_dofs, p_var*n_u_dofs, n_u_dofs, n_p_dofs);

      Kvu.reposition (v_var*n_v_dofs, u_var*n_v_dofs, n_v_dofs, n_u_dofs);
      Kvv.reposition (v_var*n_v_dofs, v_var*n_v_dofs, n_v_dofs, n_v_dofs);
      Kvp.reposition (v_var*n_v_dofs, p_var*n_v_dofs, n_v_dofs, n_p_dofs);

      Kpu.reposition (p_var*n_u_dofs, u_var*n_u_dofs, n_p_dofs, n_u_dofs);
      Kpv.reposition (p_var*n_u_dofs, v_var*n_u_dofs, n_p_dofs, n_v_dofs);
      Kpp.reposition (p_var*n_u_dofs, p_var*n_u_dofs, n_p_dofs, n_p_dofs);

      Fu.reposition (u_var*n_u_dofs, n_u_dofs);
      Fv.reposition (v_var*n_u_dofs, n_v_dofs);
      Fp.reposition (p_var*n_u_dofs, n_p_dofs);

      // Now we will build the element matrix.
      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
        {
          // Assemble the u-velocity row
          // uu coupling
          for (unsigned int i=0; i<n_u_dofs; i++)
            for (unsigned int j=0; j<n_u_dofs; j++)
              Kuu(i,j) += JxW[qp]*(dphi[i][qp]*dphi[j][qp]);

          // up coupling
          for (unsigned int i=0; i<n_u_dofs; i++)
            for (unsigned int j=0; j<n_p_dofs; j++)
              Kup(i,j) += -JxW[qp]*psi[j][qp]*dphi[i][qp](0);


          // Assemble the v-velocity row
          // vv coupling
          for (unsigned int i=0; i<n_v_dofs; i++)
            for (unsigned int j=0; j<n_v_dofs; j++)
              Kvv(i,j) += JxW[qp]*(dphi[i][qp]*dphi[j][qp]);

          // vp coupling
          for (unsigned int i=0; i<n_v_dofs; i++)
            for (unsigned int j=0; j<n_p_dofs; j++)
              Kvp(i,j) += -JxW[qp]*psi[j][qp]*dphi[i][qp](1);


          // Assemble the pressure row
          // pu coupling
          for (unsigned int i=0; i<n_p_dofs; i++)
            for (unsigned int j=0; j<n_u_dofs; j++)
              Kpu(i,j) += -JxW[qp]*psi[i][qp]*dphi[j][qp](0);

          // pv coupling
          for (unsigned int i=0; i<n_p_dofs; i++)
            for (unsigned int j=0; j<n_v_dofs; j++)
              Kpv(i,j) += -JxW[qp]*psi[i][qp]*dphi[j][qp](1);

        } // end of the quadrature point qp-loop

      // At this point the interior element integration has
      // been completed.  However, we have not yet addressed
      // boundary conditions.  For this example we will only
      // consider simple Dirichlet boundary conditions imposed
      // via the penalty method. The penalty method used here
      // is equivalent (for Lagrange basis functions) to lumping
      // the matrix resulting from the L2 projection penalty
      // approach introduced in example 3.
      {
        // The following loops over the sides of the element.
        // If the element has no neighbor on a side then that
        // side MUST live on a boundary of the domain.
        for (auto s : elem->side_index_range())
          if (elem->neighbor_ptr(s) == libmesh_nullptr)
            {
              std::unique_ptr<const Elem> side (elem->build_side_ptr(s));

              // Loop over the nodes on the side.
              for (auto ns : side->node_index_range())
                {
                  // The location on the boundary of the current
                  // node.

                  // const Real xf = side->point(ns)(0);
                  const Real yf = side->point(ns)(1);

                  // The penalty value.  \f$ \frac{1}{\epsilon \f$
                  const Real penalty = 1.e10;

                  // The boundary values.

                  // Set u = 1 on the top boundary, 0 everywhere else
                  const Real u_value = (yf > .99) ? 1. : 0.;

                  // Set v = 0 everywhere
                  const Real v_value = 0.;

                  // Find the node on the element matching this node on
                  // the side.  That defined where in the element matrix
                  // the boundary condition will be applied.
                  for (auto n : elem->node_index_range())
                    if (elem->node_id(n) == side->node_id(ns))
                      {
                        // Matrix contribution.
                        Kuu(n,n) += penalty;
                        Kvv(n,n) += penalty;

                        // Right-hand-side contribution.
                        Fu(n) += penalty*u_value;
                        Fv(n) += penalty*v_value;
                      }
                } // end face node loop
            } // end if (elem->neighbor(side) == libmesh_nullptr)
      } // end boundary condition section

      // If this assembly program were to be used on an adaptive mesh,
      // we would have to apply any hanging node constraint equations.
      dof_map.constrain_element_matrix_and_vector (Ke, Fe, dof_indices);

      // The element matrix and right-hand-side are now built
      // for this element.  Add them to the global matrix and
      // right-hand-side vector.  The NumericMatrix::add_matrix()
      // and NumericVector::add_vector() members do this for us.
      system.matrix->add_matrix (Ke, dof_indices);
      system.rhs->add_vector    (Fe, dof_indices);
    } // end of element loop
}
