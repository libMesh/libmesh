// The libMesh Finite Element Library.
// Copyright (C) 2002-2020 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



// <h1>Reduced Basis: Example 6 - Heat transfer on a curved domain in 3D</h1>
// \author David Knezevic
// \date 2012
//
// In this example we consider heat transfer modeled by a Poisson
// equation with Robin boundary condition:
//
//  -kappa \Laplacian u = 1, on \Omega
//  -kappa du\dn = kappa Bi u, on \partial\Omega_Biot,
//  u = 0 on \partial\Omega_Dirichlet,
//
// We consider a reference domain \Omega_hat =
// [-0.2,0.2]x[-0.2,0.2]x[0,3], and the physical domain is then obtain
// via the parametrized mapping:
//
//  x = -1/mu + (1/mu+x_hat)*cos(mu*z_hat)
//  y = y_hat
//  z = (1/mu+x_hat)*sin(mu*z_hat)
//
// for (x_hat,y_hat,z_hat) \in \Omega_hat. (Here "hats" denotes
// reference domain.)  Also, the "reference Dirichlet boundaries" are
// [-0.2,0.2]x[-0.2,0.2]x{0} and [-0.2,0.2]x[-0.2,0.2]x{3}, and the
// remaining boundaries are the "Biot" boundaries.
//
// Then, after putting the PDE into weak form and mapping it to the
// reference domain, we obtain:
//  \kappa \int_\Omega_hat
//    [ (1+mu*x_hat) v_x w_x + (1+mu*x_hat) v_y w_y + 1/(1+mu*x_hat) v_z w_z ]
//    + \kappa Bi \int_\partial\Omega_hat_Biot1 (1-0.2mu) u v
//    + \kappa Bi \int_\partial\Omega_hat_Biot2 (1+mu x_hat) u v
//    + \kappa Bi \int_\partial\Omega_hat_Biot3 (1+0.2mu) u v
//    = \int_\Omega_hat (1+mu x_hat) v
// where
//  \partial\Omega_hat_Biot1 = [-0.2] x [-0.2,0.2] x [0,3]
//  \partial\Omega_hat_Biot2 = [-0.2,0.2] x {-0.2} x [0,3] \UNION [-0.2,0.2] x {0.2} x [0,3]
//  \partial\Omega_hat_Biot3 = [0.2] x [-0.2,0.2] x [0,3]
//
// The term
//  \kappa \int_\Omega_hat 1/(1+mu*x_hat) v_z w_z
// is "non-affine" (in the Reduced Basis sense), since we can't
// express it in the form \sum theta_q(kappa,mu) a(v,w). As a result,
// (as in reduced_basis_ex4) we must employ the Empirical
// Interpolation Method (EIM) in order to apply the Reduced Basis
// method here.
//
// The approach we use is to construct an EIM approximation, G_EIM, to
// the vector-valued function G(x_hat,y_hat;mu) = (1 + mu*x_hat, 1 +
// mu*x_hat, 1/(1+mu*x_hat)) and then we express the "volumetric
// integral part" of the left-hand side operator as a(v,w;mu) =
// \int_\hat\Omega G_EIM(x_hat,y_hat;mu) \dot (v_x w_x, v_y w_y, v_z
// w_z).  (We actually only need EIM for the third component of
// G_EIM, but it's helpful to demonstrate "vector-valued" EIM here.)

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
#include "eim_classes.h"
#include "assembly.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;

// Define a function to scale the mesh according to the parameter.
void transform_mesh_and_plot(EquationSystems & es,
                             Real curvature,
                             const std::string & filename);

// The main program.
int main (int argc, char ** argv)
{
  // Initialize libMesh.
  LibMeshInit init (argc, argv);

  // This example requires a linear solver package.
  libmesh_example_requires(libMesh::default_solver_package() != INVALID_SOLVER_PACKAGE,
                           "--enable-petsc, --enable-trilinos, or --enable-eigen");

#if !defined(LIBMESH_HAVE_XDR)
  // We need XDR support to write out reduced bases
  libmesh_example_requires(false, "--enable-xdr");
#elif defined(LIBMESH_DEFAULT_SINGLE_PRECISION)
  // XDR binary support requires double precision
  libmesh_example_requires(false, "--disable-singleprecision");
#elif defined(LIBMESH_DEFAULT_TRIPLE_PRECISION)
  // I have no idea why long double isn't working here... [RHS]
  libmesh_example_requires(false, "double precision");
#elif defined(LIBMESH_ENABLE_BLOCKED_STORAGE)
  // This example dies with "New nonzero at (0,2) caused a malloc"
  // when blocked storage is enabled.
  libmesh_example_requires(false, "--disable-blocked-storage");
#endif

  // This is a 3D example
  libmesh_example_requires(3 == LIBMESH_DIM, "3D support");

  // This example requires libmesh to be configured with both
  // DirichletBoundary and Cap'n Proto support.
#if !defined(LIBMESH_ENABLE_DIRICHLET) || !defined(LIBMESH_HAVE_CAPNPROTO)
  libmesh_example_requires(false, "--enable-dirichlet --enable-capnp");
#else

  // Parse the input file using GetPot
  std::string eim_parameters = "eim.in";
  std::string rb_parameters  = "rb.in";
  std::string main_parameters = "reduced_basis_ex6.in";
  GetPot infile(main_parameters);

  unsigned int n_elem_xy = infile("n_elem_xy", 1);
  unsigned int n_elem_z  = infile("n_elem_z", 1);

  // Do we write the RB basis functions to disk?
  bool store_basis_functions = infile("store_basis_functions", true);

  // Read the "online_mode" flag from the command line
  GetPot command_line (argc, argv);
  int online_mode = 0;
  if (command_line.search(1, "-online_mode"))
    online_mode = command_line.next(online_mode);

  // Create a mesh, with dimension to be overridden by build_cube, on
  // the default MPI communicator.  We currently have to create a
  // ReplicatedMesh here due to a reduced_basis regression with
  // DistributedMesh
  ReplicatedMesh mesh(init.comm());

  MeshTools::Generation::build_cube (mesh,
                                     n_elem_xy, n_elem_xy, n_elem_z,
                                     -0.2, 0.2,
                                     -0.2, 0.2,
                                     0., 3.,
                                     HEX8);

  // Create an equation systems object.
  EquationSystems equation_systems (mesh);

  SimpleEIMConstruction & eim_construction =
    equation_systems.add_system<SimpleEIMConstruction> ("EIM");
  SimpleRBConstruction & rb_construction =
    equation_systems.add_system<SimpleRBConstruction> ("RB");

  // Initialize the data structures for the equation system.
  equation_systems.init ();

  // Print out some information about the "truth" discretization
  equation_systems.print_info();
  mesh.print_info();

  // Initialize the standard RBEvaluation object
  SimpleRBEvaluation rb_eval(mesh.comm());

  // Initialize the EIM RBEvaluation object
  SimpleEIMEvaluation eim_rb_eval(mesh.comm());

  // Set the rb_eval objects for the RBConstructions
  eim_construction.set_rb_eim_evaluation(eim_rb_eval);
  rb_construction.set_rb_evaluation(rb_eval);

  if (!online_mode) // Perform the Offline stage of the RB method
    {
      // Read data from input file and print state
      eim_construction.process_parameters_file(eim_parameters);
      eim_construction.print_info();

      // Perform the EIM Greedy and write out the data
      eim_construction.initialize_eim_construction();
      eim_construction.train_eim_approximation();

#if defined(LIBMESH_HAVE_CAPNPROTO)
      RBDataSerialization::RBEIMEvaluationSerialization rb_eim_eval_writer(eim_rb_eval);
      rb_eim_eval_writer.write_to_file("rb_eim_eval.bin");
#else
      eim_construction.get_rb_evaluation().legacy_write_offline_data_to_files("eim_data");
#endif

      // Read data from input file and print state
      rb_construction.process_parameters_file(rb_parameters);

      // attach the EIM theta objects to the RBEvaluation
      eim_rb_eval.initialize_eim_theta_objects();
      rb_eval.get_rb_theta_expansion().attach_multiple_A_theta(eim_rb_eval.get_eim_theta_objects());

      // attach the EIM assembly objects to the RBConstruction
      eim_construction.initialize_eim_assembly_objects();
      rb_construction.get_rb_assembly_expansion().attach_multiple_A_assembly(eim_construction.get_eim_assembly_objects());

      // Print out the state of rb_construction now that the EIM objects have been attached
      rb_construction.print_info();

      // Need to initialize _after_ EIM greedy so that
      // the system knows how many affine terms there are
      rb_construction.initialize_rb_construction();
      rb_construction.train_reduced_basis();

#if defined(LIBMESH_HAVE_CAPNPROTO)
      RBDataSerialization::RBEvaluationSerialization rb_eval_writer(rb_construction.get_rb_evaluation());
      rb_eval_writer.write_to_file("rb_eval.bin");
#else
      rb_construction.get_rb_evaluation().legacy_write_offline_data_to_files("rb_data");
#endif

      // Write out the basis functions, if requested
      if (store_basis_functions)
        {
          // Write out the basis functions
          rb_construction.get_rb_evaluation().write_out_basis_functions(rb_construction,
                                                                        "rb_data");
        }
    }
  else // Perform the Online stage of the RB method
    {
#if defined(LIBMESH_HAVE_CAPNPROTO)
      RBDataDeserialization::RBEIMEvaluationDeserialization rb_eim_eval_reader(eim_rb_eval);
      rb_eim_eval_reader.read_from_file("rb_eim_eval.bin");
#else
      eim_rb_eval.legacy_read_offline_data_from_files("eim_data");
#endif

      // attach the EIM theta objects to rb_eval objects
      eim_rb_eval.initialize_eim_theta_objects();
      rb_eval.get_rb_theta_expansion().attach_multiple_A_theta(eim_rb_eval.get_eim_theta_objects());

      // Read in the offline data for rb_eval
#if defined(LIBMESH_HAVE_CAPNPROTO)
      RBDataDeserialization::RBEvaluationDeserialization rb_eval_reader(rb_eval);
      rb_eval_reader.read_from_file("rb_eval.bin", /*read_error_bound_data*/ true);
#else
      rb_eval.legacy_read_offline_data_from_files("rb_data");
#endif

      // Get the parameters at which we will do a reduced basis solve
      Real online_curvature = infile("online_curvature", 0.);
      Real online_Bi        = infile("online_Bi", 0.);
      Real online_kappa     = infile("online_kappa", 0.);
      RBParameters online_mu;
      online_mu.set_value("curvature", online_curvature);
      online_mu.set_value("Bi", online_Bi);
      online_mu.set_value("kappa", online_kappa);
      rb_eval.set_parameters(online_mu);
      rb_eval.print_parameters();
      rb_eval.rb_solve(rb_eval.get_n_basis_functions());

      // plot the solution, if requested
      if (store_basis_functions)
        {
          // read in the data from files
          rb_eval.read_in_basis_functions(rb_construction, "rb_data");

          rb_construction.load_rb_solution();

          transform_mesh_and_plot(equation_systems, online_curvature, "RB_sol.e");
        }
    }

#endif // LIBMESH_ENABLE_DIRICHLET

  return 0;
}

#ifdef LIBMESH_ENABLE_DIRICHLET

void transform_mesh_and_plot(EquationSystems & es,
                             Real curvature,
                             const std::string & filename)
{
  // Loop over the mesh nodes and move them!
  MeshBase & mesh = es.get_mesh();

  for (auto & node : mesh.node_ptr_range())
    {
      Real x = (*node)(0);
      Real z = (*node)(2);

      (*node)(0) = -1./curvature + (1./curvature + x)*cos(curvature*z);
      (*node)(2) = (1./curvature + x)*sin(curvature*z);
    }

#ifdef LIBMESH_HAVE_EXODUS_API
  ExodusII_IO(mesh).write_equation_systems(filename, es);
#endif
}

#endif // LIBMESH_ENABLE_DIRICHLET
