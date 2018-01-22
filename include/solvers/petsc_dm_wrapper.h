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

#ifndef LIBMESH_PETSC_DM_WRAPPER_H
#define LIBMESH_PETSC_DM_WRAPPER_H

#include "libmesh/libmesh_common.h"

#ifdef LIBMESH_HAVE_PETSC

#include <vector>
#include <memory>

// PETSc includes
#include <petsc.h>

namespace libMesh
{

/**
 * This class defines a wrapper around the PETSc DM infrastructure.
 * By coordinating DM data structures with libMesh, we can use libMesh
 * mesh hierarchies for geometric multigrid. Additionally, by setting the
 * DM data, we can additionally (with or without multigrid) define recursive
 * fieldsplits of our variables.
 *
 * \author Paul T. Bauman, Boris Boutkov
 * \date 2018
 */
class PetscDMWrapper
{
public:

  PetscDMWrapper() = default;

  ~PetscDMWrapper(){};

private:

  //! Vector of DMs for all grid levels
  std::vector<std::unique_ptr<DM>> _dms;

  //! Vector of PETScSections for all grid levels
  std::vector<std::unique_ptr<PetscSection>> _sections;

  //! Vector of star forests for all grid levels
  std::vector<std::unique_ptr<PetscSF>> _star_forests;

  //! Init all the n_mesh_level dependent data structures
  void init_dm_data(unsigned int n_levels);

  //! Get reference to DM for the given mesh level
  /**
   * init_dm_data() should be called before this function.
   */
  DM & get_dm(unsigned int level)
  { libmesh_assert(level < _dms.size());
    return *(_dms[level].get()); }

  //! Get reference to PetscSection for the given mesh level
  /**
   * init_dm_data() should be called before this function.
   */
  PetscSection & get_section(unsigned int level)
  { libmesh_assert(level < _sections.size());
    return *(_sections[level].get()); }

  //! Get reference to PetscSF for the given mesh level
  /**
   * init_dm_data() should be called before this function.
   */
  PetscSF & get_star_forest(unsigned int level)
  { libmesh_assert(level < _star_forests.size());
    return *(_star_forests[level].get()); }

};

}

#endif // #ifdef LIBMESH_HAVE_PETSC

#endif // LIBMESH_PETSC_DM_WRAPPER_H
