// $Id: enum_preconditioner_type.h,v 1.4 2004-10-14 22:02:47 jwpeterson Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2004  Benjamin S. Kirk, John W. Peterson
  
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



#ifndef __enum_preconditioner_type_h__
#define __enum_preconditioner_type_h__

// C++ includes

// Local includes
#include "libmesh_config.h"




/*
 * The \p libMeshEnums namespace is the namespace all \p enum definitions
 * should be put into.
 */

// ------------------------------------------------------------
// enum PreconditionerType definition
namespace libMeshEnums {
  
  /**
   * Defines an \p enum for preconditioner types
   */
  enum PreconditionerType {IDENTITY_PRECOND =0,
			   JACOBI_PRECOND,
			   BLOCK_JACOBI_PRECOND,
			   SOR_PRECOND,
			   SSOR_PRECOND,
			   EISENSTAT_PRECOND,
			   ASM_PRECOND,
			   CHOLESKY_PRECOND,
			   ICC_PRECOND,
			   ILU_PRECOND,
			   LU_PRECOND,
			   USER_PRECOND,
			   SHELL_PRECOND,
			   
			   INVALID_PRECONDITIONER};
}

using namespace libMeshEnums;



#endif




