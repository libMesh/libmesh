// $Id: ex1.C,v 1.3 2003-03-26 13:55:21 benkirk Exp $

// The Next Great Finite Element Library.
// Copyright (C) 2003  Benjamin S. Kirk
  
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




/**
 * C++ include files that we need
 */
#include <iostream>

/**
 * Functions to initialize the library.
 */
#include "libmesh.h"

/**
 * Basic include files needed for the mesh functionality.
 */
#include "mesh.h"



/**
 * \mainpage Example 1
 *
 * \section Introduction
 *
 * This is the first example program.  It simply demonstrates
 * how to create a mesh object.  A mesh is read from file,
 * information is printed to the screen, and the mesh is then
 * written 
 */
int main (int argc, char** argv)
{
  /**
   * Initialize the library.  This is necessary because the library
   * may depend on a number of other libraries (i.e. MPI  and Petsc)
   * that require initialization before use.
   */
  libMesh::init (argc, argv);

  /**
   * Force all our objects to have local scope.  By declaring
   * libMesh objects in the next pair of braces we can assert
   * that they will go out of scope (and should have been deleted)
   * before we return from main.  This allows the library to do
   * internal reference counting and assure memory is not leaked.
   */  
  {    
    /**
     * Check for proper usage. The program is designed to be run
     * as follows:
     *
     * ./ex1 -d DIM input_mesh_name [output_mesh_name]
     *
     * where [output_mesh_name] is an optional parameter giving
     * a filename to write the mesh into.
     */
    if (argc < 4)
      {
	std::cerr << "Usage: " << argv[0] << " -d 2 in.mesh [out.mesh]"
		  << std::endl;
	
	/**
	 * This handy function will print the file name, line number,
	 * and then abort.  Currrently the library does not use C++
	 * exception handling.
	 */
	error();
      }
    
    /**
     * Get the dimensionality of the mesh from argv[2]
     */
    const unsigned int dim = atoi(argv[2]);
    
    /**
     * Create a mesh with the requested dimension.
     */
    Mesh mesh(dim);
    
    /**
     * Read the input mesh.
     */
    mesh.read (argv[3]);
    
    /**
     * Print information about the mesh to the screen.
     */
    mesh.print_info();
    
    /**
     * Write the output mesh if the user specified an
     * output file name.
     */
    if (argc == 5)
      mesh.write (argv[4]);

    /**
     * At this closing brace all of our objects will be forced
     * out of scope and will get deconstructed.
     */
  }

  /**
   * All done.  Call the libMesh::close() function to close any
   * external libraries and check for leaked memory.  To be absolutey
   * certain this is called last we will return its value.  This
   * also allows main to return nonzero if memory is leaked, which
   * can be useful for testing purposes.
   */
  return libMesh::close();
}
