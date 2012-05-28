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

// <h1>Reduced Basis: Example 5 - Reduced Basis Cantilever</h1>
// Reduced Basis version of systems_of_equations_ex4: 2D cantilever

// In this example we consider the same problem as systems_of_equations_ex4,
// but we introduce one parameter, which is the thickness of the cantilever.
// (Note that for simplicity we do not consider a rigorous lower bound for
// the coercivity constant --- to compute this bound one can use the rbOOmit
// SCM classes.)
//
// We consider three parameters in this problem:
//  y_scaling: scale the mesh in the y-direction
//  x_load: the traction in the x-direction on the right boundary of the cantilever
//  y_load: the traction in the y-direction on the right boundary of the cantilever

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

// Define a function to scale the mesh according to the parameter.
void scale_mesh_and_plot(EquationSystems& es, const RBParameters& mu, const std::string& filename);

// The main program.
int main (int argc, char** argv)
{
  // Initialize libMesh.
  LibMeshInit init (argc, argv);

  const unsigned int dim = 2;

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
  libmesh_example_assert(dim <= LIBMESH_DIM, "2D support");

  std::string parameters_filename = "reduced_basis_ex5.in";
  GetPot infile(parameters_filename);

  unsigned int n_elem_x = infile("n_elem_x", 1);
  unsigned int n_elem_y = infile("n_elem_y", 1);
  Real x_size           = infile("x_size", 1.);
  Real y_size           = infile("y_size", 1.);

  bool store_basis_functions = infile("store_basis_functions", true); // Do we write the RB basis functions to disk?

  // Read the "online_mode" flag from the command line
  GetPot command_line (argc, argv);
  int online_mode = 0;
  if ( command_line.search(1, "-online_mode") )
    online_mode = command_line.next(online_mode);

  // Build a mesh.
  Mesh mesh (dim);
  MeshTools::Generation::build_square (mesh,
                                       n_elem_x, n_elem_y,
                                       0., x_size,
                                       0., y_size,
                                       QUAD9);

  // Create an equation systems object.
  EquationSystems equation_systems (mesh);

  // Build a new RBEvaluation object which will be used to perform
  // Reduced Basis calculations. This is required in both the
  // "Offline" and "Online" stages.
  SimpleRBEvaluation rb_eval;

  if(!online_mode) // Perform the Offline stage of the RB method
  {
    // We override RBConstruction with SimpleRBConstruction in order to
    // specialize a few functions for this particular problem.
    SimpleRBConstruction & rb_con =
      equation_systems.add_system<SimpleRBConstruction> ("RBElasticity");

    // Initialize the data structures for the equation system.
    equation_systems.init ();

    // Print out some information about the "truth" discretization
    equation_systems.print_info();
    mesh.print_info();
  
    // We need to give the RBConstruction object a pointer to
    // our RBEvaluation object
    rb_con.set_rb_evaluation(rb_eval);

    // Read in the data that defines this problem from the specified text file
    rb_con.process_parameters_file(parameters_filename);

    // Print out info that describes the current setup of rb_con
    rb_con.print_info();

    // Prepare rb_con for the Construction stage of the RB method.
    // This sets up the necessary data structures and performs
    // initial assembly of the "truth" affine expansion of the PDE.
    rb_con.initialize_rb_construction();

    // Compute the reduced basis space by computing "snapshots", i.e.
    // "truth" solves, at well-chosen parameter values and employing
    // these snapshots as basis functions.
    rb_con.train_reduced_basis();
    
    // Write out the data that will subsequently be required for the Evaluation stage
    rb_con.get_rb_evaluation().write_offline_data_to_files();
    
    // If requested, write out the RB basis functions for visualization purposes
    if(store_basis_functions)
    {
      // If we want to be able to visualize the solution in the online stage,
      // then we should also save the state of the equation_systems object
      // so we can initialize it properly in the online stage
      equation_systems.write("equation_systems.dat", WRITE);
      
      // Write out the basis functions
      rb_con.get_rb_evaluation().write_out_basis_functions(rb_con);
    }
  }
  else // Perform the Online stage of the RB method
  {
    // Read in the reduced basis data
    rb_eval.read_offline_data_from_files();
    
    // Iinitialize online parameters
    rb_eval.initialize_parameters(parameters_filename);
    rb_eval.print_parameters();

    // Now do the Online solve using the precomputed reduced basis
    rb_eval.rb_solve( rb_eval.get_n_basis_functions() );

    if(store_basis_functions)
    {
      // initialize the EquationSystems object by reading in the state that
      // was written out in the offline stage
      equation_systems.read("equation_systems.dat", READ);
      RBConstruction& rb_con = equation_systems.get_system<RBConstruction>("RBElasticity");
      rb_con.set_rb_evaluation(rb_eval);

      // Read in the basis functions
      rb_eval.read_in_basis_functions(rb_con);
      
      // Plot the solution
      rb_con.load_rb_solution();
#ifdef LIBMESH_HAVE_EXODUS_API
      const RBParameters& rb_eval_params = rb_eval.get_parameters();
      scale_mesh_and_plot(equation_systems, rb_eval_params, "RB_displacement.e");
#endif
    }
  }

  return 0;
}

void scale_mesh_and_plot(EquationSystems& es, const RBParameters& mu, const std::string& filename)
{
  // Loop over the mesh nodes and move them!
  MeshBase& mesh = es.get_mesh();

  MeshBase::const_node_iterator       node_it  = mesh.nodes_begin();
  const MeshBase::const_node_iterator node_end = mesh.nodes_end();
  
  for( ; node_it != node_end; node_it++)
  {
    Node* node = *node_it;
    
    Real y = (*node)(1);

    (*node)(1) = y*mu.get_value("y_scaling");
  }

#ifdef LIBMESH_HAVE_EXODUS_API
  ExodusII_IO (mesh).write_equation_systems (filename, es);
#endif
  
  // Loop over the mesh nodes and move them!
  node_it = mesh.nodes_begin();
  
  for( ; node_it != node_end; node_it++)
  {
    Node* node = *node_it;

    (*node)(1) = 1./mu.get_value("y_scaling");
  }
}
