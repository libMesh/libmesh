// $Id: enum_solver_package.h,v 1.3 2003-02-13 22:56:07 benkirk Exp $

// The Next Great Finite Element Library.
// Copyright (C) 2002  Benjamin S. Kirk, John W. Peterson
  
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

// C++ includes

// Local includes
#include "mesh_config.h"




/*
 * The \p MeshEnums namespace is the namespace all \p enum definitions
 * should be put into.
 */

// ------------------------------------------------------------
// enum SolverType definition
namespace MeshEnums {
  
  /**
   * Defines an \p enum for various linear solver packages.
   * This allows for run-time switching between solver packages
   * 
   */
  enum SolverPackage
    { 
#ifdef HAVE_PETSC
      PETSC_SOLVERS,
#endif
#ifdef HAVE_LASPACK
      LASPACK_SOLVERS,
#endif
      
      INVALID_SOLVER_PACKAGE
    };
}

using namespace MeshEnums;



#endif




