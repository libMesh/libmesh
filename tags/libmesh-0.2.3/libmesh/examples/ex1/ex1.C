// $Id: ex1.C,v 1.1 2003-01-30 19:15:23 benkirk Exp $

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
   * Check for proper usage.
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
    };

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
   * All done.  
   */
  return 0;
};
