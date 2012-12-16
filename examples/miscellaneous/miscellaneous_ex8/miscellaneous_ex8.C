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



// <h1>Miscellaneous Example 7 - Can use the PetscDMNonlinearSolver (available in PETSc-3.3.0 or above) to solve a VI version of the problem.</h1>
//
// LibMesh interfaces directly with PETSc's variational inequality solver,
// this example shows how to do it.

// Example include files
#include "libmesh/libmesh.h"
#include "libmesh/meshfree_interpolation.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;



int main(int argc, char** argv)
{
  // Initialize libMesh.
  LibMeshInit init (argc, argv);
  {
    InverseDistanceInterpolation idi;

    std::cout << idi;
  }
  return 0;
}
