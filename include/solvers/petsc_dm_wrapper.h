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
  // Forward declarations
  class System;
  class DofMap;

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

  PetscDMWrapper(){};

  ~PetscDMWrapper(){};

  void init_and_attach_petscdm(const System & system, SNES & snes);

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

  //! Takes System, empty PetscSection and populates the PetscSection
  /**
   * Take the System in its current take and an empty PetscSection and then
   * populate the PetscSection. The PetscSection is comprised of global "point"
   * numbers, where a "point" in PetscDM parlance is a geometric entity: node, edge,
   * face, or element. Then, we also add the DoF numbering for each variable
   * for each of the "points". The PetscSection, together the with PetscSF
   * will allow for recursive fieldsplits from the command line using PETSc.
   *
   * TODO: Currently only populates nodal DoF data. Need to generalize
   *       to edge, face, and interior DoFs as well.
   */
  void build_section(const System & system, PetscSection & section);

  //! Takes System, empty PetscSF and populates the PetscSF
  /**
   * The PetscSF (star forest) is a cousin of PetscSection. PetscSection
   * has the DoF info, and PetscSF gives the parallel distribution of the
   * DoF info. So PetscSF should only be necessary when we have more than
   * one MPI rank. Essentially, we are copying the DofMap.send_list(): we
   * are specifying the local dofs, what rank communicates that dof info
   * (for off-processor dofs that are communicated) and the dofs local
   * index on that rank.
   *
   * https://jedbrown.org/files/StarForest.pdf
   */
  void build_sf( const System & system, PetscSF & star_forest );

  //! Helper function for build_section.
  /**
   * This function will set the id for each "point" on the current processor within
   * the PetscSection. Should be O(n_active_local_elems).
   */
  void set_point_range_in_section( const System & system, PetscSection & section);

  //! Helper function for build_section.
  /**
   * This function will set the DoF info for each "point" in the PetscSection.
   * Should be O(n_active_local_elems).
   */
  void add_dofs_to_section( const System & system, PetscSection & section );

  //! Helper function to sanity check PetscSection construction
  /**
   * I needed this originally to help check parallel cases and the (what I think is)
   * the redundant call to PetscSectionSetDof within add_dofs_to_section.
   */
  void check_section_n_dofs_match( const System & system, PetscSection & section );

  //! Find the MPI rank of the degree of freedom, dof.
  /**
   * Helper function for build_sf
   */
  PetscInt find_dof_rank( unsigned int dof, const DofMap & dof_map ) const;

};

}

#endif // #ifdef LIBMESH_HAVE_PETSC

#endif // LIBMESH_PETSC_DM_WRAPPER_H
