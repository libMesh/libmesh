/* $Id: ex2.C,v 1.15 2004-11-15 22:09:06 benkirk Exp $ */

/* The Next Great Finite Element Library. */
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

 // <h1>Example 2 - Defining a Simple System</h1>
 //
 // This is the second example program.  It demonstrates how to
 // create an equation system for a simple scalar system.  This
 // example will also introduce some of the issues involved with using Petsc
 // in your application.
 //
 // This is the first example program that indirectly
 // uses the Petsc library.  By default equation data is stored
 // in Petsc vectors, which may span multiple processors.  Before
 // Petsc is used it must be initialized via libMesh::init().  Note that
 // by passing argc and argv to Petsc you may specify
 // command line arguments to Petsc.  For example, you might
 // try running this example as:
 //
 // ./ex2 -log_info
 //
 // to see what Petsc is doing behind the scenes or
 //
 // ./ex2 -log_summary
 // 
 // to get a summary of what Petsc did.
 // Among other things, libMesh::init() initializes the MPI
 // communications library on your system if you haven't already
 // done so.

// C++ include files that we need
#include <iostream>
//Basic include file needed for the mesh functionality.
#include "libmesh.h"
#include "mesh.h"
// Include file that defines various mesh generation utilities
#include "mesh_generation.h"
// Include file that defines (possibly multiple) systems of equations.
#include "equation_systems.h"
// Include file that defines a simple steady system
#include "transient_system.h"

int main (int argc, char** argv)
{
  libMesh::init (argc, argv);
  // This set of braces are used to force object scope.  This
  // way we can guarantee all our objects are destroyed before calling
  // libMesh::close() at the end of the program.
  {
    // A brief message to the user to inform her of the
    // exact name of the program being run, and its command line.
    std::cout << "Running " << argv[0];
    for (int i=1; i<argc; i++)
      std::cout << " " << argv[i];
    std::cout << std::endl << std::endl;
    
    // Create a 2D mesh.
    Mesh mesh (2);
    
    // Use the MeshTools::Generation mesh generator to create a uniform
    // grid on the unit square.  By default a mesh of QUAD4
    // elements will be created.  We instruct the mesh generator
    // to build a mesh of 5x5 elements.
    MeshTools::Generation::build_cube (mesh, 5, 5);
    
    // Print information about the mesh to the screen.
    mesh.print_info();

    // Create an equation systems object. This object can
    // contain multiple systems of different 
    // flavors for solving loosely coupled physics.  Each system can 
    // contain multiple variables of different approximation orders.  
    // Here we will simply create a single system with one variable.  
    // Later on, other flavors of systems will be introduced.  For the 
    // moment, we use the general system.
    // The EquationSystems object needs a reference to the mesh
    // object, so the order of construction here is important.
    EquationSystems equation_systems (mesh);
    
    // Add a flag "test" that is visible for all systems.  This
    // helps in inter-system communication.
    equation_systems.set_flag ("test");
      
    // Set a simulation-specific parameter visible for all systems.
    // This helps in inter-system-communication.
    equation_systems.set_parameter ("dummy") = 42.;
      
    // Set another simulation-specific parameter 
    equation_systems.set_parameter ("nobody") = 0.;
    
    // Now we declare the system and its variables.
    // We begin by adding a "SteadyStytem" to the
    // EquationSystems object, and we give it the name
    // "Simple System".
    equation_systems.add_system<TransientImplicitSystem> ("Simple System");
      
    // Adds the variable "u" to "Simple System".  "u"
    // will be approximated using first-order approximation.
    equation_systems("Simple System").add_variable("u", FIRST);
      
    // Initialize the data structures for the equation system.
    equation_systems.init();
      
    // Prints information about the system to the screen.
    equation_systems.print_info();

    // Write the equation system if the user specified an
    // output file name.  Note that there are two possible
    // formats to write to.  Specifying libMeshEnums::WRITE will
    // create a formatted ASCII file.  Optionally, you can specify
    // libMeshEnums::ENCODE and get an XDR-encoded binary file.
    //
    // We will write the data, clear the object, and read the file
    // we just wrote.  This is simply to demonstrate capability.
    // Note that you might use this in an application to periodically
    // dump the state of your simulation.  You can then restart from
    // this data later.
    if (argc == 2)
      if (argv[1][0] != '-')
	{
	  std::cout << "<<< Writing system to file " << argv[1]
		    << std::endl;
	  
	  // Write the system.
	  equation_systems.write (argv[1], libMeshEnums::WRITE);
	  
	  // Clear the equation systems data structure.
	  equation_systems.clear ();

	  std::cout << ">>> Reading system from file " << argv[1]
		    << std::endl << std::endl;
	  
	  // Read the file we just wrote.  This better
	  // work!
	  equation_systems.read (argv[1], libMeshEnums::READ);

	  // Print the information again.
	  equation_systems.print_info();
	}
    
    // All our objects will be destroyed at this closing brace.
    // That way we can safely call PetscFinalize() and be sure we
    // don't have any Petsc-dependent objects lurking around!
  }

  // Call libMesh::close(), which in turn finalizes PETSc and/or
  // MPI.
  return libMesh::close();
}
