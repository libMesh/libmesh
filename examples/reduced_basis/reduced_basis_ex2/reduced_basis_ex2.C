// The libMesh Finite Element Library.
// Copyright (C) 2002-2019 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

// rbOOmit: An implementation of the Certified Reduced Basis method.
// Copyright (C) 2009, 2010 David J. Knezevic
// This file is part of rbOOmit.



// <h1>Reduced Basis Example 2 - Successive Constraint Method</h1>
// \author David Knezevic
// \date 2011
//
// In this example we extend reduced_basis_ex1 to solve a steady
// convection-diffusion problem on the unit square via the Reduced
// Basis Method. In this case, we modify the PDE so that it no longer
// has a parameter-independent coercivity constant. Therefore, in
// order to obtain an error bound, we need to employ the Successive
// Constraint Method (SCM) implemented in
// RBSCMConstruction/RBSCMEvaluation to obtain a parameter-dependent
// lower bound for the coercivity constant.
//
// The PDE being solved is div(k*grad(u)) + Beta*grad(u) = f
// k is the diffusion coefficient :
// - constant in the domain 0<=x<0.5 , its value is given by the first parameter mu[0]
// - constant in the domain 0.5<=x<=1 , its value is given by the second parameter mu[1]
// Beta is the convection velocity :
// - constant in the whole domain
// - equal to zero in the y-direction
// - its value in the x-direction is given by the third (and last) parameter mu[2]
// Boundary conditions :
// - dyu=0 on top and bottom
// - u=0 on the left side
// - dxu + Beta*u = 0 on the right side

// C++ include files that we need
#include <iostream>
#include <algorithm>
#include <cstdlib> // *must* precede <cmath> for proper std:abs() on PGI, Sun Studio CC
#include <cmath>
#include <set>

// Basic include file needed for the mesh functionality.
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/equation_systems.h"
#include "libmesh/dof_map.h"
#include "libmesh/getpot.h"
#include "libmesh/elem.h"
#include "libmesh/rb_data_serialization.h"
#include "libmesh/rb_data_deserialization.h"
#include "libmesh/enum_solver_package.h"

// local includes
#include "rb_classes.h"
#include "assembly.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;

// The main program.
int main (int argc, char ** argv)
{
  // Initialize libMesh.
  LibMeshInit init (argc, argv);

  // This example requires SLEPc and GLPK
#if !defined(LIBMESH_HAVE_SLEPC) || !defined(LIBMESH_HAVE_GLPK)
  libmesh_example_requires(false, "--enable-slepc --enable-glpk");
#else

#if !defined(LIBMESH_HAVE_XDR)
  // We need XDR support to write out reduced bases
  libmesh_example_requires(false, "--enable-xdr");
#elif defined(LIBMESH_DEFAULT_SINGLE_PRECISION)
  // XDR binary support requires double precision
  libmesh_example_requires(false, "--disable-singleprecision");
#endif
  // FIXME: This example currently segfaults with Trilinos?
  libmesh_example_requires(libMesh::default_solver_package() == PETSC_SOLVERS, "--enable-petsc");

  // FIXME: This example needs lapack, which is serial only
  libmesh_example_requires(init.comm().size() == 1, "mpirun -np 1");

  // Skip this 2D example if libMesh was compiled as 1D-only.
  libmesh_example_requires(2 <= LIBMESH_DIM, "2D support");

  // Parse the input file (reduced_basis_ex2.in) using GetPot
  std::string parameters_filename = "reduced_basis_ex2.in";
  GetPot infile(parameters_filename);

  unsigned int n_elem = infile("n_elem", 1);       // Determines the number of elements in the "truth" mesh
  const unsigned int dim = 2;                      // The number of spatial dimensions

  bool store_basis_functions = infile("store_basis_functions", true); // Do we write the RB basis functions to disk?

  // Read the "online_mode" flag from the command line
  GetPot command_line (argc, argv);
  int online_mode = 0;
  if (command_line.search(1, "-online_mode"))
    online_mode = command_line.next(online_mode);

  // Build a mesh on the default MPI communicator.
  Mesh mesh (init.comm(), dim);
  MeshTools::Generation::build_square (mesh,
                                       n_elem, n_elem,
                                       0., 1.,
                                       0., 1.,
                                       QUAD4);

  // Create an equation systems object.
  EquationSystems equation_systems (mesh);

  // We override RBConstruction with SimpleRBConstruction in order to
  // specialize a few functions for this particular problem.
  SimpleRBConstruction & rb_con =
    equation_systems.add_system<SimpleRBConstruction> ("RBConvectionDiffusion");

  // Initialize the SCM Construction object
  RBSCMConstruction & rb_scm_con =
    equation_systems.add_system<RBSCMConstruction> ("RBSCMConvectionDiffusion");
  rb_scm_con.set_RB_system_name("RBConvectionDiffusion");
  rb_scm_con.add_variable("p", FIRST);

  // Initialize the data structures for the equation system.
  equation_systems.init ();

  // Print out some information about the "truth" discretization
  equation_systems.print_info();
  mesh.print_info();

  // Set parameters for the eigenvalue problems that will be solved by rb_scm_con
  equation_systems.parameters.set<unsigned int>("eigenpairs")    = 1;
  equation_systems.parameters.set<unsigned int>("basis vectors") = 3;
  equation_systems.parameters.set<unsigned int>
    ("linear solver maximum iterations") = 1000;

  // Build a new RBEvaluation object which will be used to perform
  // Reduced Basis calculations. This is required in both the
  // "Offline" and "Online" stages.
  SimpleRBEvaluation rb_eval(mesh.comm());

  // We need to give the RBConstruction object a pointer to
  // our RBEvaluation object
  rb_con.set_rb_evaluation(rb_eval);

  // We also need a SCM evaluation object to perform SCM calculations
  RBSCMEvaluation rb_scm_eval(mesh.comm());
  rb_scm_eval.set_rb_theta_expansion(rb_eval.get_rb_theta_expansion());

  // Tell rb_eval about rb_scm_eval
  rb_eval.rb_scm_eval = &rb_scm_eval;

  // Finally, need to give rb_scm_con and rb_eval a pointer to the
  // SCM evaluation object, rb_scm_eval
  rb_scm_con.set_rb_scm_evaluation(rb_scm_eval);

  if (!online_mode) // Perform the Offline stage of the RB method
    {
      // Read in the data that defines this problem from the specified text file
      rb_con.process_parameters_file(parameters_filename);
      rb_scm_con.process_parameters_file(parameters_filename);

      // Print out info that describes the current setup of rb_con
      rb_con.print_info();
      rb_scm_con.print_info();

      // Prepare rb_con for the Construction stage of the RB method.
      // This sets up the necessary data structures and performs
      // initial assembly of the "truth" affine expansion of the PDE.
      rb_con.initialize_rb_construction();

      // Perform the SCM Greedy algorithm to derive the data required
      // for rb_scm_eval to provide a coercivity lower bound.
      rb_scm_con.perform_SCM_greedy();

      // Compute the reduced basis space by computing "snapshots", i.e.
      // "truth" solves, at well-chosen parameter values and employing
      // these snapshots as basis functions.
      rb_con.train_reduced_basis();

      // Write out the data that will subsequently be required for the Evaluation stage
#if defined(LIBMESH_HAVE_CAPNPROTO)
      RBDataSerialization::RBEvaluationSerialization rb_eval_writer(rb_con.get_rb_evaluation());
      rb_eval_writer.write_to_file("rb_eval.bin");
      RBDataSerialization::RBSCMEvaluationSerialization rb_scm_eval_writer(rb_scm_con.get_rb_scm_evaluation());
      rb_scm_eval_writer.write_to_file("rb_scm_eval.bin");
#else
      rb_con.get_rb_evaluation().legacy_write_offline_data_to_files("rb_data");
      rb_scm_con.get_rb_scm_evaluation().legacy_write_offline_data_to_files("scm_data");
#endif

      // If requested, write out the RB basis functions for visualization purposes
      if (store_basis_functions)
        {
          // Write out the basis functions
          rb_con.get_rb_evaluation().write_out_basis_functions(rb_con, "rb_data");
        }
    }
  else // Perform the Online stage of the RB method
    {

      // Read in the reduced basis data
#if defined(LIBMESH_HAVE_CAPNPROTO)
      RBDataDeserialization::RBEvaluationDeserialization rb_eval_reader(rb_eval);
      rb_eval_reader.read_from_file("rb_eval.bin", /*read_error_bound_data*/ true);
      RBDataDeserialization::RBSCMEvaluationDeserialization rb_scm_eval_reader(rb_scm_eval);
      rb_scm_eval_reader.read_from_file("rb_scm_eval.bin");
#else
      rb_eval.legacy_read_offline_data_from_files("rb_data");
      rb_scm_eval.legacy_read_offline_data_from_files("scm_data");
#endif

      // Read in online_N and initialize online parameters
      unsigned int online_N = infile("online_N", 1);
      Real online_mu_0 = infile("online_mu_0", 0.);
      Real online_mu_1 = infile("online_mu_1", 0.);
      Real online_mu_2 = infile("online_mu_2", 0.);
      RBParameters online_mu;
      online_mu.set_value("mu_0", online_mu_0);
      online_mu.set_value("mu_1", online_mu_1);
      online_mu.set_value("mu_2", online_mu_2);
      rb_eval.set_parameters(online_mu);
      rb_eval.print_parameters();

      // Now do the Online solve using the precomputed reduced basis
      rb_eval.rb_solve(online_N);

      // Print out outputs as well as the corresponding output error bounds.
      libMesh::out << "output 1, value = " << rb_eval.RB_outputs[0]
                   << ", bound = " << rb_eval.RB_output_error_bounds[0]
                   << std::endl;
      libMesh::out << "output 2, value = " << rb_eval.RB_outputs[1]
                   << ", bound = " << rb_eval.RB_output_error_bounds[1]
                   << std::endl;
      libMesh::out << "output 3, value = " << rb_eval.RB_outputs[2]
                   << ", bound = " << rb_eval.RB_output_error_bounds[2]
                   << std::endl;
      libMesh::out << "output 4, value = " << rb_eval.RB_outputs[3]
                   << ", bound = " << rb_eval.RB_output_error_bounds[3]
                   << std::endl << std::endl;

      if (store_basis_functions)
        {
          // Read in the basis functions
          rb_eval.read_in_basis_functions(rb_con, "rb_data");

          // Plot the solution
          rb_con.load_rb_solution();
#ifdef LIBMESH_HAVE_EXODUS_API
          ExodusII_IO(mesh).write_equation_systems ("RB_sol.e", equation_systems);
#endif

          // Plot the first basis function that was generated from the train_reduced_basis
          // call in the Offline stage
          rb_con.load_basis_function(0);
#ifdef LIBMESH_HAVE_EXODUS_API
          ExodusII_IO(mesh).write_equation_systems ("bf0.e", equation_systems);
#endif
        }
    }

  return 0;

#endif // LIBMESH_HAVE_SLEPC && LIBMESH_HAVE_GLPK
}
