/* rbOOmit: An implementation of the Certified Reduced Basis method. */
/* Copyright (C) 2009, 2010 David J. Knezevic */
/*     This file is part of rbOOmit. */

/* rbOOmit is free software; you can redistribute it and/or */
/* modify it under the terms of the GNU Lesser General Public */
/* License as published by the Free Software Foundation; either */
/* version 2.1 of the License, or (at your option) any later version. */
  
/* rbOOmit is distributed in the hope that it will be useful, */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU */
/* Lesser General Public License for more details. */
  
/* You should have received a copy of the GNU Lesser General Public */
/* License along with this library; if not, write to the Free Software */
/* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA */

// <h1>Reduced Basis: Example 6 - Poisson equation on a curved domain in 3D</h1>

// In this example we consider the Poisson equation -\Laplacian u = 1 with
// zero Dirichlet boundary conditions on a domain with parametrized curvature.
// We consider a reference domain [0,0.2] x [0,0.2] x [0,1], and apply the parametrized
// transformation:
//  F_0(x_hat,y_hat,z_hat;mu) = -1/mu + (1/mu+x_hat)*cos(mu*z_hat)
//  F_1(x_hat,y_hat,z_hat;mu) = y_hat
//  F_2(x_hat,y_hat,z_hat;mu) = (1/mu+x_hat)*sin(mu*z_hat)
// to obtain the "physical domain" on which we solve the PDE. Here, x = F(x_hat),
// and the "hat" denotes reference domain.

// This geometric mapping yields a "non-affine" parametric PDE, hence (as in
// reduced_basis_ex4) we must employ the Empirical Interpolation Method (EIM)
// in order to apply the Reduced Basis method.

// Note that it turns out that det(J_F) = 1 + mu * x_hat, hence the right-hand
// side functional is given by:
//  f(v;mu) = \int_\hat\Omega v + mu \int_\hat\Omega \hat x v
// which is in fact "affine" (in the RB sense). As a result, we only require EIM
// to construct the approximation of the left-hand side operator.

// For the left-hand side operator, we have:
//  a(v,w;mu) = \int_\hat\Omega [ (1+mu*x_hat) v_x w_x + v_y w_y + 1/(1+mu*x_hat) v_z w_z ]
// which is non-affine due to the 1/(1+mu*x_hat) term. The approach we use is to
// construct an EIM approximation, G_EIM, to the vector-valued function
//  G(x_hat,y_hat;mu) = (1 + mu*x_hat, 1/(1+mu*x_hat))
// and then
//  a(v,w;mu) = \int_\hat\Omega G_EIM(x_hat,y_hat;mu) \dot (v_x w_x, v_z w_z) + \int_\hat\Omega v_y w_y
// (We could use EIM to only approximate the second component of G, but we want to demonstrate
// a "vector-valued EIM" case in this example.)

// C++ include files that we need
#include <iostream>
#include <algorithm>
#include <cmath>
#include <set>

// Basic include file needed for the mesh functionality.
#include "libmesh.h"
#include "mesh.h"
#include "mesh_generation.h"
#include "exodusII_io.h"
#include "equation_systems.h"
#include "dof_map.h"
#include "getpot.h"
#include "elem.h"

// local includes
#include "rb_classes.h"
#include "eim_classes.h"
#include "assembly.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;

// Define a function to scale the mesh according to the parameter.
void transform_mesh_and_plot(EquationSystems& es, Real curvature, const std::string& filename);

// The main program.
int main (int argc, char** argv)
{
  // Initialize libMesh.
  LibMeshInit init (argc, argv);

#if !defined(LIBMESH_HAVE_XDR)
  // We need XDR support to write out reduced bases
  libmesh_example_assert(false, "--enable-xdr");
#elif defined(LIBMESH_DEFAULT_SINGLE_PRECISION)
  // XDR binary support requires double precision
  libmesh_example_assert(false, "--disable-singleprecision");
#endif
  // FIXME: This example currently segfaults with Trilinos?
  libmesh_example_assert(libMesh::default_solver_package() == PETSC_SOLVERS, "--enable-petsc");

  // This is a 3D example
  libmesh_example_assert(3 == LIBMESH_DIM, "3D support");
  
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
  if ( command_line.search(1, "-online_mode") )
    online_mode = command_line.next(online_mode);

  // Build a mesh.
  Mesh mesh;
  MeshTools::Generation::build_cube (mesh,
                                     n_elem_xy, n_elem_xy, n_elem_z,
                                     0., 0.2,
                                     0., 0.2,
                                     0., 1.,
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
  SimpleRBEvaluation rb_eval;

  // Initialize the EIM RBEvaluation object
  SimpleEIMEvaluation eim_rb_eval;
  
  // Set the rb_eval objects for the RBConstructions
  eim_construction.set_rb_evaluation(eim_rb_eval);
  rb_construction.set_rb_evaluation(rb_eval);

  if(!online_mode) // Perform the Offline stage of the RB method
  {
    // Read data from input file and print state
    eim_construction.process_parameters_file(eim_parameters);
    eim_construction.print_info();
  
    // Perform the EIM Greedy and write out the data
    eim_construction.initialize_rb_construction();
    
    eim_construction.train_reduced_basis();
    eim_construction.get_rb_evaluation().write_offline_data_to_files("eim_data");
    
    // Read data from input file and print state
    rb_construction.process_parameters_file(rb_parameters);

    // attach the EIM theta objects to the RBConstruction and RBEvaluation objects
    eim_rb_eval.initialize_eim_theta_objects();
    rb_eval.get_rb_theta_expansion().attach_multiple_A_theta(eim_rb_eval.get_eim_theta_objects());
    
    // attach the EIM assembly objects to the RBConstruction object
    eim_construction.initialize_eim_assembly_objects();
    rb_construction.get_rb_assembly_expansion().attach_multiple_A_assembly(eim_construction.get_eim_assembly_objects());

    // Print out the state of rb_construction now that the EIM objects have been attached
    rb_construction.print_info();

    // Need to initialize _after_ EIM greedy so that
    // the system knows how many affine terms there are
    rb_construction.initialize_rb_construction();
    rb_construction.train_reduced_basis();
    rb_construction.get_rb_evaluation().write_offline_data_to_files("rb_data");

    // Write out the basis functions, if requested
    if(store_basis_functions)
    {
      // Write out the basis functions
      eim_construction.get_rb_evaluation().write_out_basis_functions(eim_construction,"eim_data");
      rb_construction.get_rb_evaluation().write_out_basis_functions(rb_construction,"rb_data");
    }
  }
  else // Perform the Online stage of the RB method
  {
    eim_rb_eval.read_offline_data_from_files("eim_data");

    // attach the EIM theta objects to rb_eval objects
    eim_rb_eval.initialize_eim_theta_objects();
    rb_eval.get_rb_theta_expansion().attach_multiple_A_theta(eim_rb_eval.get_eim_theta_objects());
    
    // Read in the offline data for rb_eval
    rb_eval.read_offline_data_from_files("rb_data");

    // Get the parameters at which we will do a reduced basis solve
    Real online_curvature = infile("online_curvature", 0.);
    RBParameters online_mu;
    online_mu.set_value("curvature", online_curvature);
    rb_eval.set_parameters(online_mu);
    rb_eval.print_parameters();
    rb_eval.rb_solve( rb_eval.get_n_basis_functions() );

    // plot the solution, if requested
    if(store_basis_functions)
    {
      // read in the data from files
      eim_rb_eval.read_in_basis_functions(eim_construction,"eim_data");
      rb_eval.read_in_basis_functions(rb_construction,"rb_data");

      eim_construction.load_rb_solution();
      rb_construction.load_rb_solution();

      transform_mesh_and_plot(equation_systems,online_curvature,"RB_sol.e");
    }
  }

  return 0;
}

void transform_mesh_and_plot(EquationSystems& es, Real curvature, const std::string& filename)
{
  // Loop over the mesh nodes and move them!
  MeshBase& mesh = es.get_mesh();

  MeshBase::node_iterator       node_it  = mesh.nodes_begin();
  const MeshBase::node_iterator node_end = mesh.nodes_end();
  
  for( ; node_it != node_end; node_it++)
  {
    Node* node = *node_it;
    
    Real x = (*node)(0);
    Real z = (*node)(2);

    (*node)(0) = -1./curvature + (1./curvature + x)*cos(curvature*z);
    (*node)(2) = (1./curvature + x)*sin(curvature*z);
  }

#ifdef LIBMESH_HAVE_EXODUS_API
  ExodusII_IO(mesh).write_equation_systems(filename, es);
#endif
}

