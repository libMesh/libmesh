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



#ifndef LIBMESH_PETSC_AUTO_FIELDSPLIT_H
#define LIBMESH_PETSC_AUTO_FIELDSPLIT_H

#include "libmesh/libmesh_config.h"

#ifdef LIBMESH_HAVE_PETSC

#include "libmesh/petsc_macro.h"

// Petsc include files.

// PCFieldSplitSetIs seems to have appeared late in the PETSc 3.1
// releases; we'll support it in 3.2 onward so we don't have to worry
// about compilation errors
#include <petscksp.h>

namespace libMesh
{
// Forward declarations
class System;

void petsc_auto_fieldsplit (PC my_pc, const System & sys);

} // namespace libMesh


#endif // #ifdef LIBMESH_HAVE_PETSC
#endif // LIBMESH_PETSC_AUTO_FIELDSPLIT_H
