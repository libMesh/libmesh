// The libMesh Finite Element Library.
// Copyright (C) 2002-2016 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_ENUM_SUBSET_SOLVE_MODE_H
#define LIBMESH_ENUM_SUBSET_SOLVE_MODE_H

// ------------------------------------------------------------
// enum SubsetSolveMode definition
namespace libMesh {

/**
 * \enum SubsetSolveMode defines an \p enum for the
 * question what happens to the dofs outside the given subset when a
 * system is solved on a subset.
 */
enum SubsetSolveMode {
  SUBSET_ZERO = 0, //!< Set dofs outside the subset to zero.
  SUBSET_COPY_RHS, //!< Set dofs outside the subset to the value of the corresponding dofs of the right hand side.
  SUBSET_DONT_TOUCH //!< Leaves dofs outside the subset unchanged.  This is fastest, but also most confusing because it abandons the property that the solution vector is (theoretically) independent of the initial guess.
};

}

#endif // LIBMESH_ENUM_SUBSET_SOLVE_MODE_H
