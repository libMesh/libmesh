// The libMesh Finite Element Library.
// Copyright (C) 2002-2020 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#include "libmesh/petsc_macro.h"

#ifdef LIBMESH_HAVE_PETSC
#if !PETSC_VERSION_LESS_THAN(3,7,3)
#if defined(LIBMESH_ENABLE_AMR) && defined(LIBMESH_HAVE_METAPHYSICL)

// libMesh includes
#include "libmesh/petsc_macro.h"
#include "libmesh/petsc_matrix.h"
#include "libmesh/petsc_vector.h"
#include "libmesh/wrapped_petsc.h"

// PETSc includes
#ifdef I
# define LIBMESH_SAW_I
#endif
#include <petsc.h>
#ifndef LIBMESH_SAW_I
# undef I // Avoid complex.h contamination
#endif

// C++ includes
#include <vector>
#include <memory>
#include <unordered_map>
#include <map>

namespace libMesh
{
  // Forward declarations
  class System;
  class DofObject;

  //! Struct to house data regarding where in the mesh hierarchy we are located.
  struct PetscDMContext
  {
    int n_dofs;
    int mesh_dim;
    DM * coarser_dm;
    DM * finer_dm;
    DM * global_dm;
    PetscMatrix<libMesh::Number> * K_interp_ptr;
    PetscMatrix<libMesh::Number> * K_sub_interp_ptr;
    PetscMatrix<libMesh::Number> * K_restrict_ptr;
    PetscVector<libMesh::Number> * current_vec;

    //! Stores local dofs for each var for use in subprojection matrixes
    std::vector<std::vector<numeric_index_type>> dof_vec;

    //! Stores subfield ids for use in subprojection matrixes on coarser DMs
    std::vector<PetscInt> subfields;

  PetscDMContext() :
    n_dofs(-12345),
      mesh_dim(-12345),
      coarser_dm(nullptr),
      finer_dm(nullptr),
      global_dm(nullptr),
      K_interp_ptr(nullptr),
      K_sub_interp_ptr(nullptr),
      K_restrict_ptr(nullptr),
      current_vec(nullptr)
    {}

  };

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

    ~PetscDMWrapper();

    //! Destroys and clears all build DM-related data
    void clear();

    void init_and_attach_petscdm(System & system, SNES & snes);

  private:

    /**
     * Vector of DMs for all grid levels. These are PETSc objects
     * created by calling DMShellCreate(), so we are responsible for
     * cleaning them up.
     */
    std::vector<WrappedPetsc<DM>> _dms;

    /**
     * Vector of PETScSections for all grid levels.  These are PETSc
     * objects which are attached to the DM by calling DMSetLocalSection().
     * While the DM takes care of destroying existing PetscSections in
     * calls to DMSetLocalSection(), it does not appear to clean up
     * PetscSections it holds when it is destroyed, so we manage their
     * lifetimes using WrappedPetsc objects.
     */
    std::vector<WrappedPetsc<PetscSection>> _sections;

    /**
     * Vector of star forests for all grid levels. These are PETSc
     * objects which are attached to the DM by calling
     * DMSetSectionSF(). The DM seems to take care of cleaning these
     * up itself as far as I can tell, so we do not try to manage
     * their lifetime in any way.
     */
    std::vector<PetscSF> _star_forests;

    /**
     * Vector of projection matrixes for all grid levels.  These are
     * C++ objects, they are cleaned up automatically by their
     * destructors.
     */
    std::vector<std::unique_ptr<PetscMatrix<Number>>> _pmtx_vec;

    /**
     * Vector of sub projection matrixes for all grid levels for
     * fieldsplit.  These are C++ objects, they are cleaned up
     * automatically by their destructors.
     */
    std::vector<std::unique_ptr<PetscMatrix<Number>>> _subpmtx_vec;

    /**
     * Vector of internal PetscDM context structs for all grid levels
     * Pointers to these C++ objects are passed to DMShellSetContext(),
     * they are cleaned up automatically by their destructors.
     */
    std::vector<PetscDMContext> _ctx_vec;

    /**
     * Vector of solution vectors for all grid levels.  These are C++
     * objects, they are cleaned up automatically by their
     * destructors.
     */
    std::vector<std::unique_ptr<PetscVector<Number>>> _vec_vec;

    //! Stores n_dofs for each grid level, to be used for projection matrix sizing
    std::vector<unsigned int> _mesh_dof_sizes;

    //! Stores n_local_dofs for each grid level, to be used for projection vector sizing
    std::vector<unsigned int> _mesh_dof_loc_sizes;

    //! Init all the n_mesh_level dependent data structures
    void init_dm_data(unsigned int n_levels, const Parallel::Communicator & comm);

    /**
     * Get reference to DM for the given mesh level.
     * init_dm_data() should be called before this function.
     */
    DM & get_dm(unsigned int level)
      { libmesh_assert_less(level, _dms.size());
        return *_dms[level]; }

    /**
     * Get reference to PetscSection for the given mesh level.
     * init_dm_data() should be called before this function.
     */
    PetscSection & get_section(unsigned int level)
      { libmesh_assert_less(level, _sections.size());
        return *_sections[level]; }

    /**
     * Get reference to PetscSF for the given mesh level.
     * init_dm_data() should be called before this function.
     */
    PetscSF & get_star_forest(unsigned int level)
      { libmesh_assert_less(level, _star_forests.size());
        return _star_forests[level]; }

    /**
     * Takes System, empty PetscSection and populates the PetscSection.
     * Take the System in its current state and an empty PetscSection and then
     * populate the PetscSection. The PetscSection is comprised of global "point"
     * numbers, where a "point" in PetscDM parlance is a geometric entity: node, edge,
     * face, or element. Then, we also add the DoF numbering for each variable
     * for each of the "points". The PetscSection, together the with PetscSF
     * will allow for recursive fieldsplits from the command line using PETSc.
     */
    void build_section(const System & system, PetscSection & section);

    /**
     * Takes System, empty PetscSF and populates the PetscSF.
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

    /**
     * Helper function for build_section.
     * This function will count how many "points" on the current processor have
     * DoFs associated with them and give that count to PETSc. We need to cache
     * a mapping between the global node id and our local count that we do in this
     * function because we will need the local number again in the add_dofs_to_section
     * function.
     */
    void set_point_range_in_section( const System & system,
                                     PetscSection & section,
                                     std::unordered_map<dof_id_type,dof_id_type> & node_map,
                                     std::unordered_map<dof_id_type,dof_id_type> & elem_map,
                                     std::map<dof_id_type,unsigned int> & scalar_map);

    /**
     * Helper function for build_section.
     * This function will set the DoF info for each "point" in the PetscSection.
     */
    void add_dofs_to_section (const System & system,
                              PetscSection & section,
                              const std::unordered_map<dof_id_type,dof_id_type> & node_map,
                              const std::unordered_map<dof_id_type,dof_id_type> & elem_map,
                              const std::map<dof_id_type,unsigned int> & scalar_map);

    /**
     * Helper function to sanity check PetscSection construction
     * The PetscSection contains local dof information. This helper function just facilitates
     * sanity checking that in fact it only has n_local_dofs.
     */
    dof_id_type check_section_n_dofs( PetscSection & section );

    // Helper function to reduce code duplication when setting dofs in section
    void add_dofs_helper (const System & system,
                          const DofObject & dof_object,
                          dof_id_type local_id,
                          PetscSection & section);
  };

}

#endif // #if LIBMESH_ENABLE_AMR && LIBMESH_HAVE_METAPHYSICL
#endif // #if PETSC_VERSION
#endif // #ifdef LIBMESH_HAVE_PETSC

#endif // LIBMESH_PETSC_DM_WRAPPER_H
