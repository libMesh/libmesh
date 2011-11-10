// $Id$

// The libMesh Finite Element Library.
// Copyright (C) 2002-2008 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner
  
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



#ifndef __enum_solver_package_h__
#define __enum_solver_package_h__

/*
 * The \p libMeshEnums namespace is the namespace all \p enum definitions
 * should be put into.
 */

// ------------------------------------------------------------
// enum SolverType definition
namespace libMeshEnums {
  
  /**
   * Defines an \p enum for various linear solver packages.
   * This allows for run-time switching between solver packages
   * 
   */
  enum SolverPackage
    { 
      PETSC_SOLVERS=0,
      TRILINOS_SOLVERS,
      LASPACK_SOLVERS,
      SLEPC_SOLVERS,
      
      INVALID_SOLVER_PACKAGE
    };
}

using namespace libMeshEnums;


#endif // #define __enum_solver_package_h__
