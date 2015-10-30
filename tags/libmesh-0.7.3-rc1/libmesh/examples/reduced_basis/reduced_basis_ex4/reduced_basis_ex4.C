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

// <h1>Reduced Basis: Example 4</h1>

// In this example problem we develop a reduced basis approximation for a parametrized
// PDE that has "non-affine" parameter dependence. This requires the use of the
// Empirical Interpolation Method (EIM).
//
// We first use EIM to construct an affine approximation to the non-affine term,
// which is a parametrized function that is a Gaussian with "center" defined
// by the two parameters (mu_1,mu_2) \in [-1,1]^2. We then employ this EIM
// approximation in order to generate a reduced basis approximation for the
// parametrized PDE: -0.05 * Laplacian(u) = f(mu_1,mu_2), with zero Dirichlet
// boundary conditions.

// Basic include file needed for the mesh functionality.
#include "libmesh.h"
#include "mesh.h"
#include "mesh_generation.h"
#include "equation_systems.h"
#include "exodusII_io.h"
#include "getpot.h"

#include "eim_classes.h"
#include "rb_classes.h"


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

  // Define the names of the input files we will read the problem properties from
  std::string eim_parameters = "eim.in";
  std::string rb_parameters  = "rb.in";
  std::string main_parameters = "reduced_basis_ex4.in";
  GetPot infile(main_parameters);

  unsigned int n_elem = infile("n_elem", 1);       // Determines the number of elements in the "truth" mesh
  const unsigned int dim = 2;                      // The number of spatial dimensions
  bool store_basis_functions = infile("store_basis_functions", false); // Do we write out basis functions?

  // Read the "online_mode" flag from the command line
  GetPot command_line (argc, argv);
  int online_mode = 0;
  if ( command_line.search(1, "-online_mode") )
    online_mode = command_line.next(online_mode);

  // Create a mesh (just a simple square)
  Mesh mesh (dim);
  MeshTools::Generation::build_square (mesh,
                                       n_elem, n_elem,
                                       -1., 1.,
                                       -1., 1.,
                                       QUAD4);

  // Initialize the EquationSystems object for this mesh and attach
  // the EIM and RB Construction objects
  EquationSystems equation_systems (mesh);

  // Initialize the standard RBEvaluation object
  SimpleRBEvaluation rb_eval;

  // Initialize the EIM RBEvaluation object
  SimpleEIMEvaluation eim_rb_eval;

  if(!online_mode)
  {
    SimpleEIMConstruction & eim_construction =
      equation_systems.add_system<SimpleEIMConstruction> ("EIM");
    SimpleRBConstruction & rb_construction =
      equation_systems.add_system<SimpleRBConstruction> ("RB");
    
    // Initialize the data structures for the equation system.
    equation_systems.init ();
  
    // Print out some information about the "truth" discretization
    mesh.print_info();
    equation_systems.print_info();

    // Read data from input file and print state
    eim_construction.process_parameters_file(eim_parameters);
    eim_construction.print_info();
  
    // Read data from input file and print state
    rb_construction.process_parameters_file(rb_parameters);
    
    // Set the rb_eval objects for the RBConstructions
    eim_construction.rb_eval = &eim_rb_eval;
    rb_construction.rb_eval = &rb_eval;
  

    // Perform the EIM Greedy and write out the data
    eim_construction.initialize_rb_construction();
    eim_construction.train_reduced_basis();
    eim_construction.rb_eval->write_offline_data_to_files("eim_data");

    // attach the EIM theta objects to the RBConstruction and RBEvaluation objects
    eim_rb_eval.initialize_rb_theta_objects();
    rb_construction.rb_theta_expansion->attach_multiple_theta_q_f(eim_rb_eval.rb_eim_theta_vector);
    rb_construction.rb_eval->rb_theta_expansion->attach_multiple_theta_q_f(eim_rb_eval.rb_eim_theta_vector);
    
    // attach the EIM assembly objects to the RBConstruction object
    eim_construction.initialize_EIM_F_objects();
    rb_construction.get_rb_assembly_expansion().attach_multiple_F_q_assembly(eim_construction.rb_eim_f_vector);

    // Print out the state of rb_construction now that the EIM objects have been attached
    rb_construction.print_info();

    // Need to initialize _after_ EIM greedy so that
    // the system knows how many affine terms there are
    rb_construction.initialize_rb_construction();
    rb_construction.train_reduced_basis();
    rb_construction.rb_eval->write_offline_data_to_files("rb_data");

    // Write out the basis functions, if requested
    if(store_basis_functions)
    {
      // If we want to be able to visualize the solution in the online stage,
      // then we should also save the state of the equation_systems object
      // so we can initialize it properly in the online stage
      equation_systems.write("equation_systems.xda", WRITE);

      // Write out the basis functions
      eim_construction.rb_eval->write_out_basis_functions(eim_construction,"eim_data");
      rb_construction.rb_eval->write_out_basis_functions(rb_construction,"rb_data");
    }
  }
  else
  {
    eim_rb_eval.read_offline_data_from_files("eim_data");

    // attach the EIM theta objects to rb_eval objects
    eim_rb_eval.initialize_rb_theta_objects();
    rb_eval.rb_theta_expansion->attach_multiple_theta_q_f(eim_rb_eval.rb_eim_theta_vector);
    
    // Read in the offline data for rb_eval
    rb_eval.read_offline_data_from_files("rb_data");

    // Get the parameters at which we will do a reduced basis solve
    unsigned int online_N = infile("online_N",1);
    unsigned int n_parameters = infile("n_parameters",1);
    std::vector<Real> online_mu_vector(n_parameters);
    for(unsigned int i=0; i<n_parameters; i++)
    {
      online_mu_vector[i] = infile("online_mu", online_mu_vector[i], i);
    }
    
    // compute the RB solution
    rb_eval.set_current_parameters(online_mu_vector);
    rb_eval.print_current_parameters();
    rb_eval.rb_solve(online_N);

    // plot the solution, if requested
    if(store_basis_functions)
    {
      // initialize the EquationSystems object by reading in the state that
      // was written out in the offline stage
      equation_systems.read("equation_systems.xda", READ);
      RBConstruction& rb_construction =
        equation_systems.get_system<RBConstruction>("RB");
      RBConstruction& eim_construction =
        equation_systems.get_system<RBConstruction>("EIM");
      rb_construction.rb_eval = &rb_eval;
      eim_construction.rb_eval = &eim_rb_eval;

      // read in the data from files
      eim_rb_eval.read_in_basis_functions(eim_construction,"eim_data");
      rb_eval.read_in_basis_functions(rb_construction,"rb_data");

      // load the EIM approximation for online_mu
      eim_rb_eval.set_current_parameters(online_mu_vector);
      eim_rb_eval.rb_solve(eim_rb_eval.get_n_basis_functions());

      eim_construction.load_rb_solution();
      rb_construction.load_rb_solution();
#ifdef LIBMESH_HAVE_EXODUS_API
      ExodusII_IO(mesh).write_equation_systems("RB_sol.e",equation_systems);
#endif
    }
  }

}
