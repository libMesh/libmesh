// The libMesh Finite Element Library.
// Copyright (C) 2002-2018 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#include "libmesh/libmesh_common.h"

#ifdef LIBMESH_HAVE_PETSC

#include <petscsf.h>

#include "libmesh/petsc_dm_wrapper.h"
#include "libmesh/auto_ptr.h" // libmesh_make_unique

namespace libMesh
{

void PetscDMWrapper::init_dm_data(unsigned int n_levels)
{
  _dms.resize(n_levels);
  _sections.resize(n_levels);
  _star_forests.resize(n_levels);

  for( unsigned int i = 0; i < n_levels; i++ )
    {
      _dms[i] = libmesh_make_unique<DM>();
      _sections[i] = libmesh_make_unique<PetscSection>();
      _star_forests[i] = libmesh_make_unique<PetscSF>();
    }
}

} // end namespace libMesh

#endif // LIBMESH_HAVE_PETSC
