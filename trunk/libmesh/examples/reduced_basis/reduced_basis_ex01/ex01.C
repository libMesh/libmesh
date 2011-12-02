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

// local includes
#include "rb_classes.h"
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
// specified in ex01.in.

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
  
  // Parse the input file (ex01.in) using GetPot
  std::string parameters_filename = "ex01.in";
  GetPot infile(parameters_filename);

  unsigned int n_elem = infile("n_elem", 1);       // Determines the number of elements in the "truth" mesh
  const unsigned int dim = 2;                      // The number of spatial dimensions
  
  bool store_basis_functions = infile("store_basis_functions", true); // Do we write the RB basis functions to disk?

  // Read the "online_mode" flag from the command line
  GetPot command_line (argc, argv);
  int online_mode = 0;
  if ( command_line.search(1, "-online_mode") )
    online_mode = command_line.next(online_mode);

  // Build a mesh.
  Mesh mesh (dim);
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


  // Initialize the data structures for the equation system.
  equation_systems.init ();

  // Print out some information about the "truth" discretization
  equation_systems.print_info();
  mesh.print_info();

  // Build a new RBEvaluation object which will be used to perform
  // Reduced Basis calculations. This is required in both the
  // "Offline" and "Online" stages.
  SimpleRBEvaluation rb_eval;
  
  // We need to give the RBConstruction object a pointer to
  // our RBEvaluation object
  rb_con.rb_eval = &rb_eval;

  // Read in the data that defines this problem from the specified text file
  rb_con.process_parameters_file(parameters_filename);

  // Print out info that describes the current setup of rb_con
  rb_con.print_info();



  if(!online_mode) // Perform the Offline stage of the RB method
  {
    // Prepare rb_con for the Construction stage of the RB method.
    // This sets up the necessary data structures and performs
    // initial assembly of the "truth" affine expansion of the PDE.
    rb_con.initialize_rb_construction();

    // Compute the reduced basis space by computing "snapshots", i.e.
    // "truth" solves, at well-chosen parameter values and employing
    // these snapshots as basis functions.
    rb_con.train_reduced_basis();
    
    // Write out the data that will subsequently be required for the Evaluation stage
    rb_con.rb_eval->write_offline_data_to_files();
    
    // If requested, write out the RB basis functions for visualization purposes
    if(store_basis_functions)
    {
      rb_con.rb_eval->write_out_basis_functions(rb_con);
    }
  }
  else // Perform the Online stage of the RB method
  {
    // Read in the reduced basis data
    rb_eval.read_offline_data_from_files();
    
    // Get the parameters at which we do a reduced basis solve
    unsigned int online_N = infile("online_N",1);
    std::vector<Real> online_mu_vector(rb_con.get_n_params());
    for(unsigned int i=0; i<rb_con.get_n_params(); i++)
    {
      online_mu_vector[i] = infile("online_mu", online_mu_vector[i], i);
    }

    // Set the parameters to online_mu_vector
    rb_eval.set_current_parameters(online_mu_vector);
    rb_eval.print_current_parameters();

    // Now do the Online solve using the precomputed reduced basis
    rb_eval.rb_solve(online_N);

    // Print out outputs as well as the corresponding output error bounds.
    std::cout << "output 1, value = " << rb_eval.RB_outputs[0]
              << ", bound = " << rb_eval.RB_output_error_bounds[0]
              << std::endl;
    std::cout << "output 2, value = " << rb_eval.RB_outputs[1]
              << ", bound = " << rb_eval.RB_output_error_bounds[1]
              << std::endl;
    std::cout << "output 3, value = " << rb_eval.RB_outputs[2]
              << ", bound = " << rb_eval.RB_output_error_bounds[2]
              << std::endl;
    std::cout << "output 4, value = " << rb_eval.RB_outputs[3]
              << ", bound = " << rb_eval.RB_output_error_bounds[3]
              << std::endl << std::endl;

    if(store_basis_functions)
    {
      // Read in the basis functions
      rb_eval.read_in_basis_functions(rb_con);
      
      // Plot the solution
      rb_con.load_rb_solution();
      ExodusII_IO(mesh).write_equation_systems ("RB_sol.e",equation_systems);
      
      // Plot the first basis function that was generated from the train_reduced_basis
      // call in the Offline stage
      rb_con.load_basis_function(0);
      ExodusII_IO(mesh).write_equation_systems ("bf0.e",equation_systems);
    }
  }

  return 0;
}

