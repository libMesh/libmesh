// $Id$
 
// The libMesh Finite Element Library.
// Copyright (C) 2002-2007  Benjamin S. Kirk, John W. Peterson

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


#include <iostream>
#include <libmesh.h>



// Prototype for generic unit test.
void unit_test ();


int main (int argc, char** argv)
{
  // Initialize the library.  This is necessary because the library
  // may depend on a number of other libraries (i.e. MPI  and Petsc)
  // that require initialization before use. 
  libMesh::init(argc, argv);

  {
    std::cout << "Running";
    for (int i=0; i<argc; i++)
      std::cout << " " << argv[i];
    std::cout << std::endl;
    
    // Run the test.  All unit tests are implemented
    // in a function called unit_test().
    // We link the tests in one at a time.
    unit_test();
  }

  return libMesh::close();

}
