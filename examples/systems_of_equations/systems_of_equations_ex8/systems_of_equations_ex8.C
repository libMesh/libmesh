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


// <h1> Systems Example 8 - Linear elasticity with contact. </h1>
//
// In this example, we consider a linear elastic model with contact. This is a nonlinear
// problem due to the contact condition.
//
// This example borrows from systems_of_equations_ex6 (linear elastic
// constitutive law) and systems_of_equations_ex7 (NonlinearImplicitSystem to
// handle nonlinearities). Note that this example could be extended to large
// deformation elasticity by using the constitutive law from systems_of_equations_ex7.
//
// The contact condition is imposed using a penalty method. Consider two contact surfaces,
// "master" and "slave". For each quadrature point qp on "master", we define a line in the normal
// direction and see where that line intersects slave. This gives us the normal distance,
// which we refer to as the "signed_distance", since if this value is negative it means
// we detected overlap between master and slave qp, and we impose a normal force (i.e.
// contact force) of magnitude penalty * signed_distance. If there is no overlap at qp,
// then the contact force is zero. We repeat this for all quadrature points on both sides
// of the contact surface.
//
// We also have to augment the sparsity pattern to handle the contributions to the Jacobian
// from the "other" side of the contact surface, and we use a similar approach to
// miscellanoues_ex9 for this.
//
// Note that the contact penalty term is non-smooth, which can cause convergence problems
// for a nonlinear solver. This is a widely recognized issue in contact simulation.

// C++ include files that we need
#include <iostream>

// libMesh includes
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
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

int main (int argc, char** argv)
{
  LibMeshInit init (argc, argv);

  // This example requires the PETSc nonlinear solvers
  libmesh_example_requires(libMesh::default_solver_package() == PETSC_SOLVERS, "--enable-petsc");

  // This example requires PETSc >= 3.5.0 since it uses
  // PetscMatrix::update_preallocation_and_zero().
#if PETSC_VERSION_LESS_THAN(3,5,0)
  libmesh_example_requires(false, "PETSc >= 3.5.0");
#endif

  GetPot infile("systems_of_equations_ex8.in");
  const std::string approx_order = infile("approx_order", "FIRST");
  const std::string fe_family = infile("fe_family", "LAGRANGE");

  const Real young_modulus = infile("Young_modulus", 1.0);
  const Real poisson_ratio = infile("poisson_ratio", 0.3);
  const Real initial_forcing_magnitude = infile("initial_forcing_magnitude", 0.001);

  const Real nonlinear_abs_tol = infile("nonlinear_abs_tol", 1.e-8);
  const Real nonlinear_rel_tol = infile("nonlinear_rel_tol", 1.e-8);
  const unsigned int nonlinear_max_its = infile("nonlinear_max_its", 50);
  const Real contact_penalty = infile("contact_penalty", 1.e2);
  const Real contact_proximity_tol = infile("contact_proximity_tol",0.1);

  const unsigned int n_solves = infile("n_solves", 10);
  const Real force_increment = infile("force_increment", 0.001);

  Mesh mesh(init.comm());
  mesh.read("systems_of_equations_ex8.exo");

  mesh.print_info();

  EquationSystems equation_systems (mesh);

  NonlinearImplicitSystem& system =
    equation_systems.add_system<NonlinearImplicitSystem> ("NonlinearElasticity");

  LinearElasticityWithContact le(system, contact_penalty, contact_proximity_tol);

  unsigned int u_var = system.add_variable("u",
                      Utility::string_to_enum<Order>   (approx_order),
                      Utility::string_to_enum<FEFamily>(fe_family));
  unsigned int v_var = system.add_variable("v",
                      Utility::string_to_enum<Order>   (approx_order),
                      Utility::string_to_enum<FEFamily>(fe_family));
  unsigned int w_var = system.add_variable("w",
                      Utility::string_to_enum<Order>   (approx_order),
                      Utility::string_to_enum<FEFamily>(fe_family));

  // Also, initialize an ExplicitSystem to store stresses
  ExplicitSystem& stress_system =
    equation_systems.add_system<ExplicitSystem> ("StressSystem");
  stress_system.add_variable("sigma_00", CONSTANT, MONOMIAL);
  stress_system.add_variable("sigma_01", CONSTANT, MONOMIAL);
  stress_system.add_variable("sigma_02", CONSTANT, MONOMIAL);
  stress_system.add_variable("sigma_11", CONSTANT, MONOMIAL);
  stress_system.add_variable("sigma_12", CONSTANT, MONOMIAL);
  stress_system.add_variable("sigma_22", CONSTANT, MONOMIAL);
  stress_system.add_variable("vonMises", CONSTANT, MONOMIAL);

  equation_systems.parameters.set<Real>         ("nonlinear solver absolute residual tolerance") = nonlinear_abs_tol;
  equation_systems.parameters.set<Real>         ("nonlinear solver relative residual tolerance") = nonlinear_rel_tol;
  equation_systems.parameters.set<unsigned int> ("nonlinear solver maximum iterations")          = nonlinear_max_its;

  system.nonlinear_solver->residual_and_jacobian_object = &le;

  equation_systems.parameters.set<Real>("young_modulus") = young_modulus;
  equation_systems.parameters.set<Real>("poisson_ratio") = poisson_ratio;
  equation_systems.parameters.set<Real>("forcing_magnitude") = initial_forcing_magnitude;

  // Attach Dirichlet boundary conditions
  std::set<boundary_id_type> clamped_boundaries;
  clamped_boundaries.insert(MIN_Z_BOUNDARY);
  clamped_boundaries.insert(MAX_Z_BOUNDARY);

  std::vector<unsigned int> uvw;
  uvw.push_back(u_var);
  uvw.push_back(v_var);
  uvw.push_back(w_var);

  ZeroFunction<Number> zero;

  system.get_dof_map().add_dirichlet_boundary
    (DirichletBoundary (clamped_boundaries, uvw, &zero));

  system.get_dof_map().attach_extra_sparsity_object(le.get_augment_sparsity());

  equation_systems.init();
  equation_systems.print_info();

  // Provide a loop here so that we can do a sequence of solves
  // where solve n gives a good starting guess for solve n+1.
  // This "continuation" approach is helpful for solving for
  // large values of "forcing_magnitude".
  for(unsigned int count=0; count<n_solves; count++)
  {
    libMesh::out << "Performing solve " << (count+1) << ", forcing_magnitude: "
      << equation_systems.parameters.get<Real>("forcing_magnitude")
      << ", contact_penalty: " << contact_penalty << std::endl;

    system.solve();
    system.nonlinear_solver->print_converged_reason();

    libMesh::out << "System solved at nonlinear iteration "
      << system.n_nonlinear_iterations()
      << " , final nonlinear residual norm: "
      << system.final_nonlinear_residual()
      << std::endl << std::endl;

    libMesh::out << "Solution l2 norm: " << system.solution->l2_norm()
      << std::endl << std::endl;

    libMesh::out << "Computing stresses..." << std::endl;

    le.compute_stresses();

    std::stringstream filename;
    filename << "solution_" << count << ".exo";
    ExodusII_IO (mesh).write_equation_systems(
      filename.str(),
      equation_systems);

    {
      // Increment the forcing magnitude
      Real previous_forcing =
        equation_systems.parameters.get<Real>("forcing_magnitude");
      equation_systems.parameters.set<Real>("forcing_magnitude") =
        previous_forcing + force_increment;
    }
  }

  return 0;
}
