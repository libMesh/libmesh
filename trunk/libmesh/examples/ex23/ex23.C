// rbOOmit: An implementation of the Certified Reduced Basis method.
// Copyright (C) 2009, 2010 David J. Knezevic
//
//     This file is part of rbOOmit.

// rbOOmit is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
  
// rbOOmit is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

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
#include "o_string_stream.h"
#include "elem.h"

#include "simple_rb.h"
#include "assembly.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;

// In this example problem we use the Certified Reduced Basis method
// to solve a steady convection-diffusion problem on the unit square.
// The reduced basis method relies on an expansion of the PDE in the
// form \sum_q=1^Q_a theta_a^q(\mu) * a^q(u,v) = \sum_q=1^Q_f theta_f^q(\mu) f^q(v)
// where theta_a, theta_f are parameter dependent functions and 
// a^q, f^q are parameter independent operators (\mu denotes a parameter).

// We first attach the parameter dependent functions and paramater
// independent operators to the RBSystem. Then in Offline mode, a
// reduced basis space is generated and written out to the directory
// "offline_data". In Online mode, the reduced basis data in "offline_data"
// is read in and used to solve the reduced problem for the parameters
// specified in ex23.in.

// We also attach four outputs to the system which are averages over certain
// subregions of the domain. In Online mode, we print out the values of these
// outputs as well as rigorous error bounds with respect to the output
// associated with the "truth" finite element discretization.

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

  // Skip this 2D example if libMesh was compiled as 1D-only.
  libmesh_example_assert(2 <= LIBMESH_DIM, "2D support");
  
  // Parse the input file (ex23.in) using GetPot
  std::string parameters_filename = "ex23.in";
  GetPot infile(parameters_filename);

  unsigned int n_elem = infile("n_elem", 1);       // Determines the number of elements in the "truth" mesh
  const unsigned int dim = 2;                      // The number of spatial dimensions
  
  bool online_mode = infile("online_mode", false);                    // Are we in Online mode?
  bool store_basis_functions = infile("store_basis_functions", true); // Do we write the RB basis functions to disk?

  // Build a mesh.
  Mesh mesh (dim);
  MeshTools::Generation::build_square (mesh,
                                       n_elem, n_elem,
                                       0., 1.,
                                       0., 1.,
                                       QUAD4);

  // Create an equation systems object.
  EquationSystems equation_systems (mesh);

  // We override RBSystem with SimpleRBSytem in order to specialize a few functions.
  SimpleRBSystem & system =
    equation_systems.add_system<SimpleRBSystem> ("RBConvectionDiffusion");


  // Initialize the data structures for the equation system.
  equation_systems.init ();

  // Print out some information about the "truth" discretization
  equation_systems.print_info();
  mesh.print_info();

  // Now that the libMesh data structures have been initialized
  // in equation_systems.init(), we can set up the Reduced Basis system.

  // Construct a default RBTheta object, evaluate() returns 1.
  RBTheta rb_theta;

  // Construct the assembly functors
  Ex23DirichletDofAssembly dirichlet_assembly;
  ThetaA0 theta_a_0; A0 A0_assembly;
  ThetaA1 theta_a_1; A1 A1_assembly;
  ThetaA2 theta_a_2; A2 A2_assembly;
  F0 F0_assembly;

  // Define a set of output assembly functions. Each output
  // evaluates the solution average over a region defined by
  // the arguments to the OutputAssembly constructor
  OutputAssembly L0(0.7,0.8,0.7,0.8);
  OutputAssembly L1(0.2,0.3,0.7,0.8);
  OutputAssembly L2(0.2,0.3,0.2,0.3);
  OutputAssembly L3(0.7,0.8,0.2,0.3);

  system.attach_dirichlet_dof_initialization(&dirichlet_assembly);
  
  // Attach the expansion of the PDE operator. The third argument
  // here refers to assembly over the boundary of the mesh, but this
  // problem only requires internal assembly and hence it is set to NULL
  system.attach_A_q(&theta_a_0, &A0_assembly);
  system.attach_A_q(&theta_a_1, &A1_assembly);
  system.attach_A_q(&theta_a_2, &A2_assembly);

  
  // Attach the expansion of the RHS
  system.attach_F_q(&rb_theta, &F0_assembly);
  
  // Attach output assembly
  system.attach_output(&rb_theta, &L0);
  system.attach_output(&rb_theta, &L1);
  system.attach_output(&rb_theta, &L2);
  system.attach_output(&rb_theta, &L3);
  
  // We reuse the operator A0 as the inner product matrix
  system.attach_inner_prod_assembly(&A0_assembly);


  // Read in the data that defines this problem from the specified text file
  system.process_parameters_file(parameters_filename);

  // Build a new RBEvaluation object which will be used to perform
  // Reduced Basis calculations. This is required in both the
  // "Offline" and "Online" stages.
  AutoPtr<RBEvaluation> rb_eval = system.build_rb_evaluation();

  if(!online_mode) // Perform the Offline stage of the RB method
  {
    // Prepare system for the Construction stage of the RB method.
    // This sets up the necessary data structures and performs
    // initial assembly of the "truth" affine expansion of the PDE.
    system.initialize_RB_system(*rb_eval);

    // Compute the reduced basis space by computing "snapshots", i.e.
    // "truth" solves, at well-chosen parameter values and employing
    // these snapshots as basis functions.
    system.train_reduced_basis();
    
    // Write out the data that will subsequently be required for the Evaluation stage
    system.rb_eval->write_offline_data_to_files();
    
    // If request, write out the RB basis functions for visualization purposes
    if(store_basis_functions)
    {
      system.rb_eval->write_out_basis_functions(system);
    }
  }
  else // Perform the Online stage of the RB method
  {
    // Use system to initialize rb_eval
    system.initialize_rb_eval_from_system(*rb_eval);

    // Read in the reduced basis data
    rb_eval->read_offline_data_from_files();
    
    // Get the parameters at which we do a reduced basis solve
    unsigned int online_N = infile("online_N",1);
    std::vector<Real> online_mu_vector(system.get_n_params());
    for(unsigned int i=0; i<system.get_n_params(); i++)
    {
      online_mu_vector[i] = infile("online_mu", online_mu_vector[i], i);
    }

    // Set the parameters to online_mu_vector
    rb_eval->set_current_parameters(online_mu_vector);
    rb_eval->print_current_parameters();

    // Now do the Online solve using the precomputed reduced basis
    rb_eval->RB_solve(online_N);

    // Print out outputs as well as the corresponding output error bounds.
    std::cout << "output 1, value = " << rb_eval->RB_outputs[0]
              << ", bound = " << rb_eval->RB_output_error_bounds[0]
              << std::endl;
    std::cout << "output 2, value = " << rb_eval->RB_outputs[1]
              << ", bound = " << rb_eval->RB_output_error_bounds[1]
              << std::endl;
    std::cout << "output 3, value = " << rb_eval->RB_outputs[2]
              << ", bound = " << rb_eval->RB_output_error_bounds[2]
              << std::endl;
    std::cout << "output 4, value = " << rb_eval->RB_outputs[3]
              << ", bound = " << rb_eval->RB_output_error_bounds[3]
              << std::endl << std::endl;

    if(store_basis_functions)
    {
      // Read in the basis functions
      rb_eval->read_in_basis_functions(system);

      // Point system to rb_eval so that we can visualize the RB solution
      // and RB basis functions on the finite element mesh.
      system.rb_eval = rb_eval.get();
      
      // Plot the solution
      system.load_RB_solution();
      ExodusII_IO(mesh).write_equation_systems ("RB_sol.e",equation_systems);
      
      // Plot the first basis function that was generated from the train_reduced_basis
      // call in the Offline stage
      system.load_basis_function(0);
      ExodusII_IO(mesh).write_equation_systems ("bf0.e",equation_systems);
    }
  }

  return 0;
}
