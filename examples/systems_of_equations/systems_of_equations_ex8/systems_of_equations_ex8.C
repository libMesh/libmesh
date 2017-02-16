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



// <h1> Systems Example 8 - "Small-sliding" contact. </h1>
// \author David Knezevic
// \date 2015
//
// In this example, we consider a linear elastic model with contact. We restrict ourselves
// to considering "small sliding", which means that we set the contact "load transfer" between
// surfaces up front, based on the undeformed mesh, and retain that load transfer throughout
// the solve. Even though we consider linear elasticity here, this is a nonlinear problem due
// to the contact condition.
//
// The contact condition is imposed using the augmented Lagrangian method, e.g. see
// Simo & Laursen (1992). For the sake of simplicity, in this example we assume that contact
// nodes are perfectly aligned (this assumption can be eliminated relatively easily).
//
// The mesh in this example consists of two disconnected cylinders. We add edge elements into
// the mesh in order to ensure correct parallel communication of data on contact surfaces,
// and also so that we do not have to manually augment the sparsity pattern.

//  We impose a displacement boundary condition of -1 in the z-direction on the top cylinder
// in order to impose contact.

// C++ include files that we need
#include <iostream>

// libMesh includes
#include "libmesh/libmesh.h"
#include "libmesh/replicated_mesh.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/equation_systems.h"
#include "libmesh/dof_map.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/getpot.h"
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/zero_function.h"
#include "libmesh/nonlinear_solver.h"
#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/petsc_macro.h"

// Local includes
#include "linear_elasticity_with_contact.h"

using namespace libMesh;

int main (int argc, char ** argv)
{
  LibMeshInit init (argc, argv);

  // This example uses an ExodusII input file
#ifndef LIBMESH_HAVE_EXODUS_API
  libmesh_example_requires(false, "--enable-exodus");
#endif

  // We use a 3D domain.
  libmesh_example_requires(LIBMESH_DIM > 2, "--disable-1D-only --disable-2D-only");

  // This example requires the PETSc nonlinear solvers
  libmesh_example_requires(libMesh::default_solver_package() == PETSC_SOLVERS, "--enable-petsc");

  GetPot infile("systems_of_equations_ex8.in");
  const std::string approx_order = infile("approx_order", "FIRST");
  const std::string fe_family = infile("fe_family", "LAGRANGE");

  const Real young_modulus = infile("Young_modulus", 1.0);
  const Real poisson_ratio = infile("poisson_ratio", 0.3);

  const Real nonlinear_abs_tol = infile("nonlinear_abs_tol", 1.e-8);
  const Real nonlinear_rel_tol = infile("nonlinear_rel_tol", 1.e-8);
  const unsigned int nonlinear_max_its = infile("nonlinear_max_its", 50);
  const Real contact_penalty = infile("contact_penalty", 1.e2);
  const Real gap_function_tol = infile("gap_function_tol", 1.e-8);

  // This example code has not been written to cope with a distributed mesh
  ReplicatedMesh mesh(init.comm());
  mesh.read("systems_of_equations_ex8.exo");

  mesh.print_info();

  EquationSystems equation_systems (mesh);

  NonlinearImplicitSystem & system =
    equation_systems.add_system<NonlinearImplicitSystem> ("NonlinearElasticity");

  LinearElasticityWithContact le(system, contact_penalty);

  unsigned int u_var =
    system.add_variable("u",
                        Utility::string_to_enum<Order>   (approx_order),
                        Utility::string_to_enum<FEFamily>(fe_family));

  unsigned int v_var =
    system.add_variable("v",
                        Utility::string_to_enum<Order>   (approx_order),
                        Utility::string_to_enum<FEFamily>(fe_family));

  unsigned int w_var =
    system.add_variable("w",
                        Utility::string_to_enum<Order>   (approx_order),
                        Utility::string_to_enum<FEFamily>(fe_family));

  // Also, initialize an ExplicitSystem to store stresses
  ExplicitSystem & stress_system =
    equation_systems.add_system<ExplicitSystem> ("StressSystem");

  stress_system.add_variable("sigma_00", CONSTANT, MONOMIAL);
  stress_system.add_variable("sigma_01", CONSTANT, MONOMIAL);
  stress_system.add_variable("sigma_02", CONSTANT, MONOMIAL);
  stress_system.add_variable("sigma_11", CONSTANT, MONOMIAL);
  stress_system.add_variable("sigma_12", CONSTANT, MONOMIAL);
  stress_system.add_variable("sigma_22", CONSTANT, MONOMIAL);
  stress_system.add_variable("vonMises", CONSTANT, MONOMIAL);

  equation_systems.parameters.set<Real> ("nonlinear solver absolute residual tolerance") = nonlinear_abs_tol;
  equation_systems.parameters.set<Real> ("nonlinear solver relative residual tolerance") = nonlinear_rel_tol;
  equation_systems.parameters.set<unsigned int> ("nonlinear solver maximum iterations") = nonlinear_max_its;

  system.nonlinear_solver->residual_and_jacobian_object = &le;

  equation_systems.parameters.set<Real>("young_modulus") = young_modulus;
  equation_systems.parameters.set<Real>("poisson_ratio") = poisson_ratio;

  // Attach Dirichlet boundary conditions
  {
    std::set<boundary_id_type> clamped_boundaries;
    clamped_boundaries.insert(MIN_Z_BOUNDARY);

    std::vector<unsigned int> uvw;
    uvw.push_back(u_var);
    uvw.push_back(v_var);
    uvw.push_back(w_var);

    ZeroFunction<Number> zero;

    // Most DirichletBoundary users will want to supply a "locally
    // indexed" functor
    system.get_dof_map().add_dirichlet_boundary
      (DirichletBoundary (clamped_boundaries, uvw, zero,
                          LOCAL_VARIABLE_ORDER));
  }
  {
    std::set<boundary_id_type> clamped_boundaries;
    clamped_boundaries.insert(MAX_Z_BOUNDARY);

    std::vector<unsigned int> uv;
    uv.push_back(u_var);
    uv.push_back(v_var);

    ZeroFunction<Number> zero;

    system.get_dof_map().add_dirichlet_boundary
      (DirichletBoundary (clamped_boundaries, uv, zero,
                          LOCAL_VARIABLE_ORDER));
  }
  {
    std::set<boundary_id_type> clamped_boundaries;
    clamped_boundaries.insert(MAX_Z_BOUNDARY);

    std::vector<unsigned int> w;
    w.push_back(w_var);

    ConstFunction<Number> neg_one(-1.);

    system.get_dof_map().add_dirichlet_boundary
      (DirichletBoundary (clamped_boundaries, w, neg_one,
                          LOCAL_VARIABLE_ORDER));
  }

  le.initialize_contact_load_paths();

  libMesh::out << "Mesh before adding edge connectors" << std::endl;
  mesh.print_info();
  le.add_contact_edge_elements();

  libMesh::out << "Mesh after adding edge connectors" << std::endl;
  mesh.print_info();

  equation_systems.init();
  equation_systems.print_info();

  libMesh::out << "Contact penalty: " << contact_penalty  << std::endl << std::endl;

  Real current_max_gap_function = std::numeric_limits<Real>::max();

  const unsigned int max_outer_iterations = 10;
  for (unsigned int outer_iteration = 0;
       outer_iteration != max_outer_iterations; ++outer_iteration)
    {
      if (current_max_gap_function <= gap_function_tol)
        {
          libMesh::out << "Outer loop converged" << std::endl;
          break;
        }

      libMesh::out << "Starting outer iteration " << outer_iteration << std::endl;

      // Perform inner iteration (i.e. Newton's method loop)
      system.solve();
      system.nonlinear_solver->print_converged_reason();

      // Perform augmented Lagrangian update
      le.update_lambdas();

      std::pair<Real, Real> least_and_max_gap_function = le.get_least_and_max_gap_function();
      Real least_gap_fn = least_and_max_gap_function.first;
      Real max_gap_fn = least_and_max_gap_function.second;

      libMesh::out << "Finished outer iteration, least gap function: "
                   << least_gap_fn
                   << ", max gap function: "
                   << max_gap_fn
                   << std::endl
                   << std::endl;

      current_max_gap_function = std::max(std::abs(least_gap_fn), std::abs(max_gap_fn));

      outer_iteration++;
    }

  libMesh::out << "Computing stresses..." << std::endl;

  le.compute_stresses();

  std::stringstream filename;
  filename << "solution.exo";
  ExodusII_IO (mesh).write_equation_systems(filename.str(),
                                            equation_systems);

  return 0;
}
