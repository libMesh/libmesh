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
// 3D cantilever beam using the Reduced Basis Method
// David Knezevic, Kyung Hoon Lee
//
// We consider four parameters in this problem:
//  x_scaling: scales the length of the cantilever
//  load_Fx: the traction in the x-direction on the right boundary of the cantilever
//  load_Fy: the traction in the y-direction on the right boundary of the cantilever
//  load_Fz: the traction in the z-direction on the right boundary of the cantilever

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
#include "assembly.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;

// Define a function to scale the mesh according to the parameter.
void scale_mesh_and_plot(EquationSystems& es, const RBParameters& mu, const std::string& filename);

// The main program.
int main(int argc, char** argv) {
	// Initialize libMesh.
	LibMeshInit init (argc, argv);

#if !defined(LIBMESH_HAVE_XDR)
  // We need XDR support to write out reduced bases
  libmesh_example_assert(false, "--enable-xdr");
#elif defined(LIBMESH_DEFAULT_SINGLE_PRECISION)
  // XDR binary support requires double precision
  libmesh_example_assert(false, "--disable-singleprecision");
#endif

  // This example only works if libMesh was compiled for 3D
  const unsigned int dim = 3;
  libmesh_example_assert(dim == LIBMESH_DIM, "3D support");

	std::string parameters_filename = "reduced_basis_ex5.in";
	GetPot infile(parameters_filename);

  unsigned int n_elem_x  = infile("n_elem_x",0);
  unsigned int n_elem_y  = infile("n_elem_y",0);
  unsigned int n_elem_z  = infile("n_elem_z",0);
  Real x_size            = infile("x_size", 0.);
  Real y_size            = infile("y_size", 0.);
  Real z_size            = infile("z_size", 0.);

	bool store_basis_functions = infile("store_basis_functions", true);

	// Read the "online_mode" flag from the command line
	GetPot command_line(argc, argv);
	int online_mode = 0;
	if ( command_line.search(1, "-online_mode") ) {
		online_mode = command_line.next(online_mode);		
	}


  Mesh mesh (dim);
  MeshTools::Generation::build_cube (mesh,
                                     n_elem_x,
                                     n_elem_y,
                                     n_elem_z,
                                     0., x_size,
                                     0., y_size,
                                     0., z_size,
                                     HEX8);
	 mesh.print_info();

	// Create an equation systems object.
	EquationSystems equation_systems(mesh);

	// We override RBConstruction with ElasticityRBConstruction in order to
	// specialize a few functions for this particular problem.
	ElasticityRBConstruction& rb_con =
		equation_systems.add_system<ElasticityRBConstruction>("RBElasticity");

	// Initialize the data structures for the equation system.
	equation_systems.init ();
  equation_systems.print_info();

	// Build a new RBEvaluation object which will be used to perform
	// Reduced Basis calculations. This is required in both the
	// "Offline" and "Online" stages.
	ElasticityRBEvaluation rb_eval;

	// We need to give the RBConstruction object a pointer to
	// our RBEvaluation object
	rb_con.set_rb_evaluation(rb_eval);

	if(!online_mode) // Perform the Offline stage of the RB method
	{
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
			// Write out the basis functions
			rb_con.get_rb_evaluation().write_out_basis_functions(rb_con);
		}
	}
	else // Perform the Online stage of the RB method
	{
		// Read in the reduced basis data
		rb_eval.read_offline_data_from_files();

		// Iinitialize online parameters
		Real online_x_scaling = infile("online_x_scaling", 0.);
		Real online_load_Fx   = infile("online_load_Fx",   0.);
		Real online_load_Fy   = infile("online_load_Fy",   0.);
		Real online_load_Fz   = infile("online_load_Fz",   0.);
		RBParameters online_mu;
		online_mu.set_value("x_scaling", online_x_scaling);
		online_mu.set_value("load_Fx",   online_load_Fx);
		online_mu.set_value("load_Fy",   online_load_Fy);
		online_mu.set_value("load_Fz",   online_load_Fz);
		rb_eval.set_parameters(online_mu);
		rb_eval.print_parameters();
		
		// Now do the Online solve using the precomputed reduced basis
		rb_eval.rb_solve( rb_eval.get_n_basis_functions() );

		if(store_basis_functions)
		{
			// Read in the basis functions
			rb_eval.read_in_basis_functions(rb_con);

			// Plot the solution
			rb_con.load_rb_solution();

			const RBParameters& rb_eval_params = rb_eval.get_parameters();
			scale_mesh_and_plot(equation_systems, rb_eval_params, "RB_sol.e");
		}		
	}

	return 0;
}

void scale_mesh_and_plot(EquationSystems& es, const RBParameters& mu, const std::string& filename)
{
  const Real x_scaling = mu.get_value("x_scaling");
  
  // Loop over the mesh nodes and move them!
  MeshBase& mesh = es.get_mesh();

  MeshBase::node_iterator       node_it  = mesh.nodes_begin();
  const MeshBase::node_iterator node_end = mesh.nodes_end();

  for( ; node_it != node_end; node_it++)
  {
    Node* node = *node_it;

    (*node)(0) *= mu.get_value("x_scaling");
  }

#ifdef LIBMESH_HAVE_EXODUS_API
  ExodusII_IO (mesh).write_equation_systems (filename, es);
#endif

  // Loop over the mesh nodes and move them!
  node_it = mesh.nodes_begin();

  for( ; node_it != node_end; node_it++)
  {
    Node* node = *node_it;

    (*node)(0) /= mu.get_value("x_scaling");
  }
}
