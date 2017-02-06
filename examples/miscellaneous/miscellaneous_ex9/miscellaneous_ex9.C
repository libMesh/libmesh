// The libMesh Finite Element Library.
// Copyright (C) 2002-2017 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



// <h1>Miscellaneous Example 9 - Implement an interface term to model
// a thermal "film resistance"</h1>
// \author David Knezevic
// \date 2013
//
// In this example we solve a Poisson problem, -\Laplacian u = f, with
// a non-standard interface condition on the domain interior which
// models a thermal "film resistance". The interface condition
// requires continuity of flux, and a jump in temperature proportional
// to the flux:
//  \nabla u_1 \cdot n = \nabla u_2 \cdot n,
//  u_1 - u_2 = R * \nabla u \cdot n
//
// To implement this PDE, we use two mesh subdomains, \Omega_1 and
// \Omega_2, with coincident boundaries, but which are not connected
// in the FE sense. Let \Gamma denote the coincident boundary.  The
// term on \Gamma takes the form:
//
//  1/R * \int_\Gamma (u_1 - u_2) (v_1 - v_2) ds,
//
// where u_1, u_2 (resp. v_1, v_2) are the trial (resp. test)
// functions on either side of \Gamma.  We implement this condition
// using C0 basis functions, but the "crack" in the mesh at \Gamma
// permits a discontinuity in the solution. We also impose a heat flux
// on the bottom surface of the mesh, and a zero Dirichlet condition
// on the top surface.
//
// In order to implement the interface condition, we need to augment
// the matrix sparsity pattern, which is handled by the class
// AugmentSparsityPatternOnInterface.


// C++ include files that we need
#include <iostream>
#include <limits>

// libMesh includes
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/mesh_refinement.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/equation_systems.h"
#include "libmesh/fe.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/dof_map.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/getpot.h"
#include "libmesh/elem.h"
#include "libmesh/fe_interface.h"
#include "libmesh/boundary_info.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/zero_function.h"
#include "libmesh/dirichlet_boundaries.h"

// local includes
#include "augment_sparsity_on_interface.h"

// define the boundary IDs in the mesh
#define MIN_Z_BOUNDARY 1
#define MAX_Z_BOUNDARY 2
#define CRACK_BOUNDARY_LOWER 3
#define CRACK_BOUNDARY_UPPER 4

// Bring in everything from the libMesh namespace
using namespace libMesh;

/**
 * Assemble the system matrix and rhs vector.
 */
void assemble_poisson(EquationSystems & es,
                      const ElementSideMap & lower_to_upper);

// The main program.
int main (int argc, char ** argv)
{
  // Initialize libMesh.
  LibMeshInit init (argc, argv);

  // This example uses an ExodusII input file
#ifndef LIBMESH_HAVE_EXODUS_API
  libmesh_example_requires(false, "--enable-exodus");
#endif

  // The sparsity augmentation code requires PETSc
  libmesh_example_requires(libMesh::default_solver_package() == PETSC_SOLVERS, "--enable-petsc");

  // Skip this 3D example if libMesh was compiled as 1D or 2D-only.
  libmesh_example_requires(3 <= LIBMESH_DIM, "3D support");

  GetPot command_line (argc, argv);

  Real R = 2.;
  if (command_line.search(1, "-R"))
    R = command_line.next(R);

  Mesh mesh(init.comm());

  EquationSystems equation_systems (mesh);

  LinearImplicitSystem & system =
    equation_systems.add_system<LinearImplicitSystem> ("Poisson");
  system.add_variable("u", FIRST, LAGRANGE);

  // We want to call assemble_poisson "manually" so that we can pass in
  // lower_to_upper, hence set assemble_before_solve = false
  system.assemble_before_solve = false;

  // Impose zero Dirichlet boundary condition on MAX_Z_BOUNDARY
  std::set<boundary_id_type> boundary_ids;
  boundary_ids.insert(MAX_Z_BOUNDARY);
  std::vector<unsigned int> variables;
  variables.push_back(0);
  ZeroFunction<> zf;

  // Most DirichletBoundary users will want to supply a "locally
  // indexed" functor
  DirichletBoundary dirichlet_bc(boundary_ids, variables, zf,
                                 LOCAL_VARIABLE_ORDER);
  system.get_dof_map().add_dirichlet_boundary(dirichlet_bc);

  // Attach an object to the DofMap that will augment the sparsity pattern
  // due to the degrees-of-freedom on the "crack"
  //
  // By attaching this object *before* reading our mesh, we also
  // ensure that the connected elements will not be deleted on a
  // distributed mesh.
  AugmentSparsityOnInterface augment_sparsity
    (mesh, CRACK_BOUNDARY_LOWER, CRACK_BOUNDARY_UPPER);
  system.get_dof_map().add_coupling_functor(augment_sparsity);

  mesh.read("miscellaneous_ex9.exo");

  equation_systems.init();
  equation_systems.print_info();

  // Set the jump term coefficient, it will be used in assemble_poisson
  equation_systems.parameters.set<Real>("R") = R;

  // Assemble and then solve
  assemble_poisson(equation_systems,
                   augment_sparsity.get_lower_to_upper());
  system.solve();

#ifdef LIBMESH_HAVE_EXODUS_API
  // Plot the solution
  ExodusII_IO (mesh).write_equation_systems ("solution.exo",
                                             equation_systems);
#endif

#ifdef LIBMESH_ENABLE_AMR
  // Possibly solve on a refined mesh next.
  MeshRefinement mesh_refinement (mesh);
  unsigned int n_refinements = 0;
  if (command_line.search("-n_refinements"))
    n_refinements = command_line.next(0);

  for (unsigned int r = 0; r != n_refinements; ++r)
    {
      std::cout << "Refining the mesh" << std::endl;

      mesh_refinement.uniformly_refine ();
      equation_systems.reinit();

      assemble_poisson(equation_systems,
                       augment_sparsity.get_lower_to_upper());
      system.solve();

#ifdef LIBMESH_HAVE_EXODUS_API
      // Plot the refined solution
      std::ostringstream out;
      out << "solution_" << r << ".exo";
      ExodusII_IO (mesh).write_equation_systems (out.str(),
                                                 equation_systems);
#endif

    }

#endif

  return 0;
}

void assemble_poisson(EquationSystems & es,
                      const ElementSideMap & lower_to_upper)
{
  const MeshBase & mesh = es.get_mesh();
  const unsigned int dim = mesh.mesh_dimension();

  Real R = es.parameters.get<Real>("R");

  LinearImplicitSystem & system = es.get_system<LinearImplicitSystem>("Poisson");

  const DofMap & dof_map = system.get_dof_map();

  FEType fe_type = dof_map.variable_type(0);

  UniquePtr<FEBase> fe (FEBase::build(dim, fe_type));
  UniquePtr<FEBase> fe_elem_face (FEBase::build(dim, fe_type));
  UniquePtr<FEBase> fe_neighbor_face (FEBase::build(dim, fe_type));

  QGauss qrule (dim, fe_type.default_quadrature_order());
  QGauss qface(dim-1, fe_type.default_quadrature_order());

  fe->attach_quadrature_rule (&qrule);
  fe_elem_face->attach_quadrature_rule (&qface);
  fe_neighbor_face->attach_quadrature_rule (&qface);

  const std::vector<Real> & JxW = fe->get_JxW();
  const std::vector<std::vector<Real> > & phi = fe->get_phi();
  const std::vector<std::vector<RealGradient> > & dphi = fe->get_dphi();

  const std::vector<Real> & JxW_face = fe_elem_face->get_JxW();

  const std::vector<Point> & qface_points = fe_elem_face->get_xyz();

  const std::vector<std::vector<Real> > & phi_face          = fe_elem_face->get_phi();
  const std::vector<std::vector<Real> > & phi_neighbor_face = fe_neighbor_face->get_phi();

  DenseMatrix<Number> Ke;
  DenseVector<Number> Fe;

  DenseMatrix<Number> Kne;
  DenseMatrix<Number> Ken;
  DenseMatrix<Number> Kee;
  DenseMatrix<Number> Knn;

  std::vector<dof_id_type> dof_indices;

  MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

  for ( ; el != end_el; ++el)
    {
      const Elem * elem = *el;

      dof_map.dof_indices (elem, dof_indices);
      const unsigned int n_dofs = dof_indices.size();

      fe->reinit (elem);

      Ke.resize (n_dofs, n_dofs);
      Fe.resize (n_dofs);

      // Assemble element interior terms for the matrix
      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
        for (unsigned int i=0; i<n_dofs; i++)
          for (unsigned int j=0; j<n_dofs; j++)
            Ke(i,j) += JxW[qp]*(dphi[i][qp]*dphi[j][qp]);

      // Boundary flux provides forcing in this example
      {
        for (unsigned int side=0; side<elem->n_sides(); side++)
          if (elem->neighbor_ptr(side) == libmesh_nullptr)
            {
              if (mesh.get_boundary_info().has_boundary_id (elem, side, MIN_Z_BOUNDARY))
                {
                  fe_elem_face->reinit(elem, side);

                  for (unsigned int qp=0; qp<qface.n_points(); qp++)
                    for (std::size_t i=0; i<phi.size(); i++)
                      Fe(i) += JxW_face[qp] * phi_face[i][qp];
                }

            }
      }

      // Add boundary terms on the crack
      {
        for (unsigned int side=0; side<elem->n_sides(); side++)
          if (elem->neighbor_ptr(side) == libmesh_nullptr)
            {
              // Found the lower side of the crack. Assemble terms due to lower and upper in here.
              if (mesh.get_boundary_info().has_boundary_id (elem, side, CRACK_BOUNDARY_LOWER))
                {
                  fe_elem_face->reinit(elem, side);

                  ElementSideMap::const_iterator ltu_it =
                    lower_to_upper.find(std::make_pair(elem, side));
                  libmesh_assert(ltu_it != lower_to_upper.end());

                  const Elem * neighbor = ltu_it->second;

                  std::vector<Point> qface_neighbor_points;
                  FEInterface::inverse_map (elem->dim(), fe->get_fe_type(),
                                            neighbor, qface_points, qface_neighbor_points);
                  fe_neighbor_face->reinit(neighbor, &qface_neighbor_points);

                  std::vector<dof_id_type> neighbor_dof_indices;
                  dof_map.dof_indices (neighbor, neighbor_dof_indices);
                  const unsigned int n_neighbor_dofs = neighbor_dof_indices.size();

                  Kne.resize (n_neighbor_dofs, n_dofs);
                  Ken.resize (n_dofs, n_neighbor_dofs);
                  Kee.resize (n_dofs, n_dofs);
                  Knn.resize (n_neighbor_dofs, n_neighbor_dofs);

                  // Lower-to-lower coupling term
                  for (unsigned int qp=0; qp<qface.n_points(); qp++)
                    for (unsigned int i=0; i<n_dofs; i++)
                      for (unsigned int j=0; j<n_dofs; j++)
                        Kee(i,j) -= JxW_face[qp] * (1./R)*(phi_face[i][qp] * phi_face[j][qp]);

                  // Lower-to-upper coupling term
                  for (unsigned int qp=0; qp<qface.n_points(); qp++)
                    for (unsigned int i=0; i<n_dofs; i++)
                      for (unsigned int j=0; j<n_neighbor_dofs; j++)
                        Ken(i,j) += JxW_face[qp] * (1./R)*(phi_face[i][qp] * phi_neighbor_face[j][qp]);

                  // Upper-to-upper coupling term
                  for (unsigned int qp=0; qp<qface.n_points(); qp++)
                    for (unsigned int i=0; i<n_neighbor_dofs; i++)
                      for (unsigned int j=0; j<n_neighbor_dofs; j++)
                        Knn(i,j) -= JxW_face[qp] * (1./R)*(phi_neighbor_face[i][qp] * phi_neighbor_face[j][qp]);

                  // Upper-to-lower coupling term
                  for (unsigned int qp=0; qp<qface.n_points(); qp++)
                    for (unsigned int i=0; i<n_neighbor_dofs; i++)
                      for (unsigned int j=0; j<n_dofs; j++)
                        Kne(i,j) += JxW_face[qp] * (1./R)*(phi_neighbor_face[i][qp] * phi_face[j][qp]);

                  system.matrix->add_matrix(Kne, neighbor_dof_indices, dof_indices);
                  system.matrix->add_matrix(Ken, dof_indices, neighbor_dof_indices);
                  system.matrix->add_matrix(Kee, dof_indices);
                  system.matrix->add_matrix(Knn, neighbor_dof_indices);
                }
            }
      }

      dof_map.constrain_element_matrix_and_vector (Ke, Fe, dof_indices);

      system.matrix->add_matrix (Ke, dof_indices);
      system.rhs->add_vector    (Fe, dof_indices);
    }
}
