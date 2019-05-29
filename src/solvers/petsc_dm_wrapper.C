// The libMesh Finite Element Library.
// Copyright (C) 2002-2019 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "libmesh/petsc_macro.h"

#ifdef LIBMESH_HAVE_PETSC
#if !PETSC_VERSION_LESS_THAN(3,7,3)
#if defined(LIBMESH_ENABLE_AMR) && defined(LIBMESH_HAVE_METAPHYSICL)

#include "libmesh/ignore_warnings.h"
#include <petscsf.h>
#if PETSC_VERSION_LESS_THAN(3,12,0)
#include <petsc/private/dmimpl.h>
#endif
#include <petscdmshell.h>
#include "libmesh/restore_warnings.h"

#include "libmesh/petsc_dm_wrapper.h"
#include "libmesh/system.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_base.h"
#include "libmesh/mesh_refinement.h"
#include "libmesh/partitioner.h"
#include "libmesh/dof_map.h"
#include "libmesh/elem.h"

namespace libMesh
{

  //--------------------------------------------------------------------
  // Functions with C linkage to pass to PETSc.  PETSc will call these
  // methods as needed.
  //
  // Since they must have C linkage they have no knowledge of a namespace.
  // We give them an obscure name to avoid namespace pollution.
  //--------------------------------------------------------------------
  extern "C"
  {

    //! Help PETSc create a subDM given a global dm when using fieldsplit
#if PETSC_VERSION_LESS_THAN(3,9,0)
    PetscErrorCode libmesh_petsc_DMCreateSubDM(DM dm, PetscInt numFields, PetscInt fields[], IS *is, DM *subdm)
#else
    PetscErrorCode libmesh_petsc_DMCreateSubDM(DM dm, PetscInt numFields, const PetscInt fields[], IS *is, DM *subdm)
#endif
    {
      PetscErrorCode ierr;
      // Basically, we copy the PETSc ShellCreateSubDM implementation,
      // but also need to set the embedding dim and also propagate
      // the relevant function pointers to the subDM for GMG purposes.
      // Since this is called by PETSc we gotta pull some of this info
      // from the context in the DM.

      // First, retrieve our context
      void * ctx = nullptr;
      ierr = DMShellGetContext(dm, & ctx);
      LIBMESH_CHKERR(ierr);
      libmesh_assert(ctx);
      PetscDMContext * p_ctx = static_cast<PetscDMContext * >(ctx);

      if (subdm)
        {
          ierr = DMShellCreate(PetscObjectComm((PetscObject) dm), subdm);
          LIBMESH_CHKERR(ierr);

          // Set the DM embedding dimension to help PetscDS (Discrete System)
          ierr = DMSetCoordinateDim(*subdm, p_ctx->mesh_dim);
          LIBMESH_CHKERR(ierr);

          // Now set the function pointers for the subDM
          // Some DMShellGet* functions only exist with PETSc >= 3.12.0.

          // Set Coarsen function pointer
#if PETSC_VERSION_LESS_THAN(3,12,0)
          if (dm->ops->coarsen)
            {
              ierr = DMShellSetCoarsen(*subdm, dm->ops->coarsen);
              LIBMESH_CHKERR(ierr);
            }
#else
          PetscErrorCode (*coarsen)(DM,MPI_Comm,DM*) = nullptr;
          ierr = DMShellGetCoarsen(dm, &coarsen);
          LIBMESH_CHKERR(ierr);
          if (coarsen)
            {
              ierr = DMShellSetCoarsen(*subdm, coarsen);
              LIBMESH_CHKERR(ierr);
            }
#endif

          // Set Refine function pointer
#if PETSC_VERSION_LESS_THAN(3,12,0)
          if (dm->ops->refine)
            {
              ierr = DMShellSetRefine(*subdm, dm->ops->refine);
              LIBMESH_CHKERR(ierr);
            }
#else
          PetscErrorCode (*refine)(DM,MPI_Comm,DM*) = nullptr;
          ierr = DMShellGetRefine(dm, &refine);
          if (refine)
            {
              ierr = DMShellSetRefine(*subdm, refine);
              LIBMESH_CHKERR(ierr);
            }
#endif

          // Set Interpolation function pointer
#if PETSC_VERSION_LESS_THAN(3,12,0)
          if (dm->ops->createinterpolation)
            {
              ierr = DMShellSetCreateInterpolation(*subdm, dm->ops->createinterpolation);
              LIBMESH_CHKERR(ierr);
            }
#else
          PetscErrorCode (*interp)(DM,DM,Mat*,Vec*) = nullptr;
          ierr = DMShellGetCreateInterpolation(*subdm, &interp);
          LIBMESH_CHKERR(ierr);
          if (interp)
            {
              ierr = DMShellSetCreateInterpolation(*subdm, interp);
              LIBMESH_CHKERR(ierr);
            }
#endif

          // Set Restriction function pointer
#if PETSC_VERSION_LESS_THAN(3,12,0)
          if (dm->ops->createrestriction)
            {
              ierr = DMShellSetCreateRestriction(*subdm, dm->ops->createrestriction);
              LIBMESH_CHKERR(ierr);
            }
#else
          PetscErrorCode (*createrestriction)(DM,DM,Mat*) = nullptr;
          ierr = DMShellGetCreateRestriction(*subdm, &createrestriction);
          LIBMESH_CHKERR(ierr);
          if (createrestriction)
            {
              ierr = DMShellSetCreateRestriction(*subdm, createrestriction);
              LIBMESH_CHKERR(ierr);
            }
#endif

          // Set CreateSubDM function pointer
#if PETSC_VERSION_LESS_THAN(3,12,0)
          if (dm->ops->createsubdm)
            {
              ierr = DMShellSetCreateSubDM(*subdm, dm->ops->createsubdm);
              LIBMESH_CHKERR(ierr);
            }
#else
          PetscErrorCode (*createsubdm)(DM,PetscInt,const PetscInt[],IS*,DM*) = nullptr;
          ierr = DMShellGetCreateSubDM(*subdm, &createsubdm);
          LIBMESH_CHKERR(ierr);
          if (createsubdm)
            {
              ierr = DMShellSetCreateSubDM(*subdm, createsubdm);
              LIBMESH_CHKERR(ierr);
            }
#endif
          // Set Context pointer
          if (ctx)
            {
              ierr = DMShellSetContext(*subdm, ctx);
              LIBMESH_CHKERR(ierr);
            }
#if PETSC_VERSION_LESS_THAN(3,11,0)
          // Lastly, Compute the subsection for the subDM
          ierr = DMCreateSubDM_Section_Private(dm, numFields, fields, is, subdm);
          LIBMESH_CHKERR(ierr);
#else
          ierr = DMCreateSectionSubDM(dm, numFields, fields, is, subdm);
          LIBMESH_CHKERR(ierr);
#endif
        }

      return 0;
    }

    //! Help PETSc identify the finer DM given a dmc
    PetscErrorCode libmesh_petsc_DMRefine(DM dmc, MPI_Comm /*comm*/, DM * dmf)
    {
      libmesh_assert(dmc);
      libmesh_assert(dmf);

      PetscErrorCode ierr;

      // extract our context from the incoming dmc
      void * ctx_c = nullptr;
      ierr = DMShellGetContext(dmc, & ctx_c);LIBMESH_CHKERR(ierr);
      libmesh_assert(ctx_c);
      PetscDMContext * p_ctx = static_cast<PetscDMContext * >(ctx_c);

      // check / set the finer DM
      libmesh_assert(p_ctx->finer_dm);
      libmesh_assert(*(p_ctx->finer_dm));
      *(dmf) = *(p_ctx->finer_dm);

      return 0;
    }

    //! Help PETSc identify the coarser DM dmc given the fine DM dmf
    PetscErrorCode libmesh_petsc_DMCoarsen(DM dmf, MPI_Comm /*comm*/, DM * dmc)
    {
      libmesh_assert(dmc);
      libmesh_assert(dmf);

      PetscErrorCode ierr;

      // Extract our context from the incoming dmf
      void * ctx_f = nullptr;
      ierr = DMShellGetContext(dmf, &ctx_f);LIBMESH_CHKERR(ierr);
      libmesh_assert(ctx_f);
      PetscDMContext * p_ctx_f = static_cast<PetscDMContext*>(ctx_f);

      // First, ensure that there exists a coarse DM that we want to
      // set. There ought to be as we created it while walking the
      // hierarchy.
      libmesh_assert(p_ctx_f->coarser_dm);
      libmesh_assert(*(p_ctx_f->coarser_dm));

      // In situations using fieldsplit we need to provide a coarser
      // DM which only has the relevant subfields in it. Since we
      // create global DMs for each mesh level, we need to also create
      // the subDMs. We do this by checking the number of fields. When
      // less than all the fields are used, we need to create the
      // proper subDMs. We get the number of fields and their names
      // from the incoming fine DM and the global reference DM
      PetscInt nfieldsf, nfieldsg;
      char ** fieldnamesf;
      char ** fieldnamesg;

      libmesh_assert(p_ctx_f->global_dm);
      DM * globaldm = p_ctx_f->global_dm;
      ierr = DMCreateFieldIS(dmf, &nfieldsf, &fieldnamesf, nullptr);
      LIBMESH_CHKERR(ierr);
      ierr = DMCreateFieldIS(*globaldm, &nfieldsg, &fieldnamesg, nullptr);
      LIBMESH_CHKERR(ierr);

      // If the probed number of fields is less than the number of
      // global fields, this amounts to PETSc 'indicating' to us we
      // are doing FS. So, we must create subDMs for the coarser
      // DMs.
      if ( nfieldsf < nfieldsg )
        {
          p_ctx_f->subfields.clear();
          p_ctx_f->subfields.resize(nfieldsf);

          // To select the subDM fields we match fine grid DM field
          //  names to their global DM counterparts. Since PETSc can
          //  internally reassign field numbering under a fieldsplit,
          //  we must extract subsections via the field names. This is
          //  admittedly gross, but c'est la vie.
          for (int i = 0; i < nfieldsf ; i++)
            {
              for (int j = 0; j < nfieldsg ;j++)
                if ( strcmp( fieldnamesg[j], fieldnamesf[i] ) == 0 )
                  p_ctx_f->subfields[i] = j;
            }

          // Next, for the found fields we create a subDM
          DM subdm;
          libmesh_petsc_DMCreateSubDM(*(p_ctx_f->coarser_dm), nfieldsf,
                                      p_ctx_f->subfields.data(), nullptr, &subdm);

          // Extract our coarse context from the created subDM so we
          // can set its subfields for use in createInterp.
          void * ctx_c = nullptr;
          ierr = DMShellGetContext(subdm, &ctx_c);
          LIBMESH_CHKERR(ierr);
          libmesh_assert(ctx_c);
          PetscDMContext * p_ctx_c = static_cast<PetscDMContext*>(ctx_c);

          // propogate subfield info to subDM
          p_ctx_c->subfields = p_ctx_f->subfields;

          // return created subDM to PETSc
          *(dmc) = subdm;
        }
      else {
        // No fieldsplit was requested so set the coarser DM to the
        // global coarser DM.
        *(dmc) = *(p_ctx_f->coarser_dm);
      }

      return 0;
    }

    //! Function to give PETSc that sets the Interpolation Matrix between two DMs
    PetscErrorCode
    libmesh_petsc_DMCreateInterpolation (DM dmc /*coarse*/, DM dmf /*fine*/,
                                         Mat * mat ,Vec * vec)
    {
      libmesh_assert(dmc);
      libmesh_assert(dmf);
      libmesh_assert(mat);
      libmesh_assert(vec); // Optional scaling (not needed for mg)

      // Get a communicator from incoming DM
      PetscErrorCode ierr;
      MPI_Comm comm;
      PetscObjectGetComm((PetscObject)dmc, &comm);

      // Extract our coarse context from the incoming DM
      void * ctx_c = nullptr;
      ierr = DMShellGetContext(dmc, &ctx_c);
      LIBMESH_CHKERR(ierr);
      libmesh_assert(ctx_c);
      PetscDMContext * p_ctx_c = static_cast<PetscDMContext*>(ctx_c);

      // Extract our fine context from the incoming DM
      void * ctx_f = nullptr;
      ierr = DMShellGetContext(dmf, &ctx_f);LIBMESH_CHKERR(ierr);
      libmesh_assert(ctx_f);
      PetscDMContext * p_ctx_f = static_cast<PetscDMContext*>(ctx_f);

      // Check for existing global projection matrix
      libmesh_assert(p_ctx_c->K_interp_ptr);

      // If were doing fieldsplit we need to construct sub projection
      // matrices. We compare the passed in number of DMs fields to a
      // global DM in order to determine if a subprojection is needed.
      PetscInt nfieldsf, nfieldsg;

      libmesh_assert(p_ctx_c->global_dm);
      DM * globaldm = p_ctx_c->global_dm;

      ierr = DMCreateFieldIS(dmf, &nfieldsf, nullptr, nullptr);
      LIBMESH_CHKERR(ierr);
      ierr = DMCreateFieldIS(*globaldm, &nfieldsg, nullptr, nullptr);
      LIBMESH_CHKERR(ierr);

      // If the probed number of fields is less than the number of
      // global fields, this amounts to PETSc 'indicating' to us we
      // are doing FS.
      if ( nfieldsf < nfieldsg)
        {
          // Loop over the fields and merge their index sets.
          std::vector<std::vector<numeric_index_type>> allrows,allcols;
          std::vector<numeric_index_type> rows,cols;
          allrows = p_ctx_f->dof_vec;
          allcols = p_ctx_c->dof_vec;

          // For internal libmesh submat extraction need to merge all
          // field dofs and then sort the vectors so that they match
          // the Projection Matrix ordering
          const int n_subfields = p_ctx_f->subfields.size();
          if ( n_subfields >= 1 )
            {
              for (int i : p_ctx_f->subfields)
                {
                  rows.insert(rows.end(), allrows[i].begin(), allrows[i].end());
                  cols.insert(cols.end(), allcols[i].begin(), allcols[i].end());
                }
              std::sort(rows.begin(),rows.end());
              std::sort(cols.begin(),cols.end());
            }

          // Now that we have merged the fine and coarse index sets
          // were ready to make the submatrix and pass it off to PETSc
          p_ctx_c->K_interp_ptr->create_submatrix (*p_ctx_c->K_sub_interp_ptr, rows, cols);

          // return to PETSc the created submatrix
          *(mat) = p_ctx_c->K_sub_interp_ptr->mat();

        } // endif less incoming DM fields than global DM fields
      else
        {
          // We are not doing fieldsplit, so return global projection
          *(mat) = p_ctx_c->K_interp_ptr->mat();
        }

      // Vec scaling isnt needed so were done.
      *(vec) = PETSC_NULL;

      return 0;
    } // end libmesh_petsc_DMCreateInterpolation

    //! Function to give PETSc that sets the Restriction Matrix between two DMs
    PetscErrorCode
    libmesh_petsc_DMCreateRestriction (DM dmc /*coarse*/, DM dmf/*fine*/, Mat * mat)
    {
      libmesh_assert(dmc);
      libmesh_assert(dmf);
      libmesh_assert(mat);

      PetscErrorCode ierr;

      // get a communicator from incoming DM
      MPI_Comm comm;
      PetscObjectGetComm((PetscObject)dmc, &comm);

      // extract our fine context from the incoming DM
      void * ctx_f = nullptr;
      ierr = DMShellGetContext(dmf, &ctx_f);LIBMESH_CHKERR(ierr);
      libmesh_assert(ctx_f);
      PetscDMContext * p_ctx_f = static_cast<PetscDMContext*>(ctx_f);

      // check / give PETSc its matrix
      libmesh_assert(p_ctx_f->K_restrict_ptr);
      *(mat) = p_ctx_f->K_restrict_ptr->mat();

      return 0;
    }

  } // end extern C functions


  PetscDMWrapper::~PetscDMWrapper()
  {
    this->clear();
  }

  void PetscDMWrapper::clear()
  {
    // PETSc will destroy the attached PetscSection, PetscSF as well as
    // other relateds such as the Projections so we just tidy up the
    // containers here.

    _dms.clear();
    _sections.clear();
    _star_forests.clear();
    _pmtx_vec.clear();
    _vec_vec.clear();
    _ctx_vec.clear();

  }

  void PetscDMWrapper::init_and_attach_petscdm(System & system, SNES & snes)
  {
    START_LOG ("init_and_attach_petscdm()", "PetscDMWrapper");

    PetscErrorCode ierr;

    MeshBase & mesh = system.get_mesh();   // Convenience
    MeshRefinement mesh_refinement(mesh); // Used for swapping between grids

    // Theres no need for these code paths while traversing the hierarchy
    mesh.allow_renumbering(false);
    mesh.allow_remote_element_removal(false);
    mesh.partitioner() = nullptr;

    // First walk over the active local elements and see how many maximum MG levels we can construct
    unsigned int n_levels = 0;
    for ( auto & elem : mesh.active_local_element_ptr_range() )
      {
        if ( elem->level() > n_levels )
          n_levels = elem->level();
      }
    // On coarse grids some processors may have no active local elements,
    // these processors shouldnt make projections
    if (n_levels >= 1)
      n_levels += 1;

    // How many MG levels did the user request?
    unsigned int usr_requested_mg_lvls = 0;
    usr_requested_mg_lvls = command_line_next("-pc_mg_levels", usr_requested_mg_lvls);

    // Only construct however many levels were requested if something was actually requested
    if ( usr_requested_mg_lvls != 0 )
      {
        // Dont request more than avail num levels on mesh, require at least 2 levels
        libmesh_assert_less_equal( usr_requested_mg_lvls, n_levels );
        libmesh_assert( usr_requested_mg_lvls > 1 );

        n_levels = usr_requested_mg_lvls;
      }
    else
      {
        // if -pc_mg_levels is not specified we just construct fieldsplit related
        // structures on the finest mesh.
        n_levels = 1;
      }


    // Init data structures: data[0] ~ coarse grid, data[n_levels-1] ~ fine grid
    this->init_dm_data(n_levels, system.comm());

    // Step 1.  contract : all active elements have no children
    mesh.contract();

    // Start on finest grid. Construct DM datas and stash some info for
    // later projection_matrix and vec sizing
    for(unsigned int level = n_levels; level >= 1; level--)
      {
        // Save the n_fine_dofs before coarsening for later projection matrix sizing
        _mesh_dof_sizes[level-1] = system.get_dof_map().n_dofs();
        _mesh_dof_loc_sizes[level-1] = system.get_dof_map().n_local_dofs();

        // Get refs to things we will fill
        DM & dm = this->get_dm(level-1);
        PetscSection & section = this->get_section(level-1);
        PetscSF & star_forest = this->get_star_forest(level-1);

        // The shell will contain other DM info
        ierr = DMShellCreate(system.comm().get(), &dm);
        LIBMESH_CHKERR(ierr);

        // Set the DM embedding dimension to help PetscDS (Discrete System)
        ierr = DMSetCoordinateDim(dm, mesh.mesh_dimension());
        CHKERRABORT(system.comm().get(),ierr);

        // Build the PetscSection and attach it to the DM
        this->build_section(system, section);
        ierr = DMSetDefaultSection(dm, section);
        LIBMESH_CHKERR(ierr);

        // We only need to build the star forest if we're in a parallel environment
        if (system.n_processors() > 1)
          {
            // Build the PetscSF and attach it to the DM
            this->build_sf(system, star_forest);
            ierr = DMSetDefaultSF(dm, star_forest);
            LIBMESH_CHKERR(ierr);
          }

        // Set PETSC's Restriction, Interpolation, Coarsen and Refine functions for the current DM
        ierr = DMShellSetCreateInterpolation ( dm, libmesh_petsc_DMCreateInterpolation );
        LIBMESH_CHKERR(ierr);

        // Not implemented. For now we rely on galerkin style restrictions
        bool supply_restriction = false;
        if (supply_restriction)
          {
            ierr = DMShellSetCreateRestriction ( dm, libmesh_petsc_DMCreateRestriction  );
            LIBMESH_CHKERR(ierr);
          }

        ierr = DMShellSetCoarsen ( dm, libmesh_petsc_DMCoarsen );
        LIBMESH_CHKERR(ierr);

        ierr = DMShellSetRefine ( dm, libmesh_petsc_DMRefine );
        LIBMESH_CHKERR(ierr);

        ierr= DMShellSetCreateSubDM(dm, libmesh_petsc_DMCreateSubDM);
        CHKERRABORT(system.comm().get(), ierr);

        // Uniformly coarsen if not the coarsest grid and distribute dof info.
        if ( level != 1 )
          {
            START_LOG ("PDM_coarsen", "PetscDMWrapper");
            mesh_refinement.uniformly_coarsen(1);
            STOP_LOG  ("PDM_coarsen", "PetscDMWrapper");

            START_LOG ("PDM_dist_dof", "PetscDMWrapper");
            system.get_dof_map().distribute_dofs(mesh);
            STOP_LOG  ("PDM_dist_dof", "PetscDMWrapper");
          }
      } // End PETSc data structure creation

    // Now fill the corresponding internal PetscDMContext for each created DM
    for( unsigned int i=1; i <= n_levels; i++ )
      {
        // Set context dimension
        (*_ctx_vec[i-1]).mesh_dim = mesh.mesh_dimension();

        // Create and attach a sized vector to the current ctx
        _vec_vec[i-1]->init( _mesh_dof_sizes[i-1] );
        _ctx_vec[i-1]->current_vec = _vec_vec[i-1].get();

        // Set a global DM to be used as reference when using fieldsplit
        _ctx_vec[i-1]->global_dm = &(this->get_dm(n_levels-1));

        if (n_levels > 1 )
          {
            // Set pointers to surrounding dm levels to help PETSc refine/coarsen
            if ( i == 1 ) // were at the coarsest mesh
              {
                (*_ctx_vec[i-1]).coarser_dm = nullptr;
                (*_ctx_vec[i-1]).finer_dm   = _dms[1].get();
              }
            else if( i == n_levels ) // were at the finest mesh
              {
                (*_ctx_vec[i-1]).coarser_dm = _dms[_dms.size() - 2].get();
                (*_ctx_vec[i-1]).finer_dm   = nullptr;
              }
            else // were in the middle of the hierarchy
              {
                (*_ctx_vec[i-1]).coarser_dm = _dms[i-2].get();
                (*_ctx_vec[i-1]).finer_dm   = _dms[i].get();
              }
          }

      } // End context creation

    // Attach a vector and context to each DM
    if ( n_levels >= 1 )
      {

        for ( unsigned int i = 1; i <= n_levels ; ++i)
          {
            DM & dm = this->get_dm(i-1);

            ierr = DMShellSetGlobalVector( dm, (*_ctx_vec[ i-1 ]).current_vec->vec() );
            LIBMESH_CHKERR(ierr);

            ierr = DMShellSetContext( dm, _ctx_vec[ i-1 ].get() );
            LIBMESH_CHKERR(ierr);
          }
      }

    // DM structures created, now we need projection matrixes if GMG is requested.
    // To prepare for projection creation go to second coarsest mesh so we can utilize
    // old_dof_indices information in the projection creation.
    if (n_levels > 1 )
      {

        // First, stash the coarse dof indices for FS purposes
        unsigned int n_vars  = system.n_vars();
        _ctx_vec[0]->dof_vec.resize(n_vars);

        for( unsigned int v = 0; v < n_vars; v++ )
          {
            std::vector<numeric_index_type> di;
            system.get_dof_map().local_variable_indices(di, system.get_mesh(), v);
            _ctx_vec[0]->dof_vec[v] = di;
          }

        START_LOG ("PDM_refine", "PetscDMWrapper");
        mesh_refinement.uniformly_refine(1);
        STOP_LOG  ("PDM_refine", "PetscDMWrapper");

        START_LOG ("PDM_dist_dof", "PetscDMWrapper");
        system.get_dof_map().distribute_dofs(mesh);
        STOP_LOG  ("PDM_dist_dof", "PetscDMWrapper");

        START_LOG ("PDM_cnstrnts", "PetscDMWrapper");
        system.reinit_constraints();
        STOP_LOG ("PDM_cnstrnts", "PetscDMWrapper");
      }

    // Create the Interpolation Matrices between adjacent mesh levels
    for ( unsigned int i = 1 ; i < n_levels ; ++i )
      {
        if ( i != n_levels )
          {
            // Stash the rest of the dof indices
            unsigned int n_vars  = system.n_vars();
            _ctx_vec[i]->dof_vec.resize(n_vars);

            for( unsigned int v = 0; v < n_vars; v++ )
              {
                std::vector<numeric_index_type> di;
                system.get_dof_map().local_variable_indices(di, system.get_mesh(), v);
                _ctx_vec[i]->dof_vec[v] = di;
              }

            unsigned int ndofs_c = _mesh_dof_sizes[i-1];
            unsigned int ndofs_f = _mesh_dof_sizes[i];

            // Create the Interpolation matrix and set its pointer
            _ctx_vec[i-1]->K_interp_ptr = _pmtx_vec[i-1].get();
            _ctx_vec[i-1]->K_sub_interp_ptr = _subpmtx_vec[i-1].get();

            unsigned int ndofs_local     = system.get_dof_map().n_dofs_on_processor(system.processor_id());
            unsigned int ndofs_old_first = system.get_dof_map().first_old_dof(system.processor_id());
            unsigned int ndofs_old_end   = system.get_dof_map().end_old_dof(system.processor_id());
            unsigned int ndofs_old_size  = ndofs_old_end - ndofs_old_first;

            // Init and zero the matrix
            _ctx_vec[i-1]->K_interp_ptr->init(ndofs_f, ndofs_c, ndofs_local, ndofs_old_size, 30 , 20);

            // Disable Mat destruction since PETSc destroys these for us
            _ctx_vec[i-1]->K_interp_ptr->set_destroy_mat_on_exit(false);
            _ctx_vec[i-1]->K_sub_interp_ptr->set_destroy_mat_on_exit(false);

            // TODO: Projection matrix sparsity pattern?
            //MatSetOption(_ctx_vec[i-1]->K_interp_ptr->mat(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);

            // Compute the interpolation matrix and set K_interp_ptr
            START_LOG ("PDM_proj_mat", "PetscDMWrapper");
            system.projection_matrix(*_ctx_vec[i-1]->K_interp_ptr);
            STOP_LOG  ("PDM_proj_mat", "PetscDMWrapper");

            // Always close matrix that contains altered data
            _ctx_vec[i-1]->K_interp_ptr->close();
          }

        // Move to next grid to make next projection
        if ( i != n_levels - 1 )
          {
            START_LOG ("PDM_refine", "PetscDMWrapper");
            mesh_refinement.uniformly_refine(1);
            STOP_LOG  ("PDM_refine", "PetscDMWrapper");

            START_LOG ("PDM_dist_dof", "PetscDMWrapper");
            system.get_dof_map().distribute_dofs(mesh);
            STOP_LOG ("PDM_dist_dof", "PetscDMWrapper");

            START_LOG ("PDM_cnstrnts", "PetscDMWrapper");
            system.reinit_constraints();
            STOP_LOG  ("PDM_cnstrnts", "PetscDMWrapper");

          }
      } // End create transfer operators. System back at the finest grid

    // Lastly, give SNES the finest level DM
    DM & dm = this->get_dm(n_levels-1);
    ierr = SNESSetDM(snes, dm);
    LIBMESH_CHKERR(ierr);

    STOP_LOG ("init_and_attach_petscdm()", "PetscDMWrapper");
  }

  void PetscDMWrapper::build_section( const System & system, PetscSection & section )
  {
    START_LOG ("build_section()", "PetscDMWrapper");

    PetscErrorCode ierr;
    ierr = PetscSectionCreate(system.comm().get(),&section);
    LIBMESH_CHKERR(ierr);

    // Tell the PetscSection about all of our System variables
    ierr = PetscSectionSetNumFields(section,system.n_vars());
    LIBMESH_CHKERR(ierr);

    // Set the actual names of all the field variables
    for( unsigned int v = 0; v < system.n_vars(); v++ )
      {
        ierr = PetscSectionSetFieldName( section, v, system.variable_name(v).c_str() );
        LIBMESH_CHKERR(ierr);
      }

    // For building the section, we need to create local-to-global map
    // of local "point" ids to the libMesh global id of that point.
    // A "point" in PETSc nomenclature is a geometric object that can have
    // dofs associated with it, e.g. Node, Edge, Face, Elem.
    // The numbering PETSc expects is continuous for the local numbering.
    // Since we're only using this interface for solvers, then we can just
    // assign whatever local id to any of the global ids. But it is local
    // so we don't need to worry about processor numbering for the local
    // point ids.
    std::unordered_map<dof_id_type,dof_id_type> node_map;
    std::unordered_map<dof_id_type,dof_id_type> elem_map;
    std::map<dof_id_type,unsigned int> scalar_map;

    // First we tell the PetscSection about all of our points that have
    // dofs associated with them.
    this->set_point_range_in_section(system, section, node_map, elem_map, scalar_map);

    // Now we can build up the dofs per "point" in the PetscSection
    this->add_dofs_to_section(system, section, node_map, elem_map, scalar_map);

    // Final setup of PetscSection
    // Until Matt Knepley finishes implementing the commented out function
    // below, the PetscSection will be assuming node-major ordering
    // so let's throw an error if the user tries to use this without
    // node-major order
    if (!libMesh::on_command_line("--node-major-dofs"))
      libmesh_error_msg("ERROR: Must use --node-major-dofs with PetscSection!");

    //else if (!system.identify_variable_groups())
    //  ierr = PetscSectionSetUseFieldOffsets(section,PETSC_TRUE);LIBMESH_CHKERR(ierr);
    //else
    //  {
    //    std::string msg = "ERROR: Only node-major or var-major ordering supported for PetscSection!\n";
    //    msg += "       var-group-major ordering not supported!\n";
    //    msg += "       Must use --node-major-dofs or set System::identify_variable_groups() = false!\n";
    //    libmesh_error_msg(msg);
    //  }

    ierr = PetscSectionSetUp(section);LIBMESH_CHKERR(ierr);

    // Sanity checking at least that local_n_dofs match
    libmesh_assert_equal_to(system.n_local_dofs(),this->check_section_n_dofs(section));

    STOP_LOG ("build_section()", "PetscDMWrapper");
  }

  void PetscDMWrapper::build_sf( const System & system, PetscSF & star_forest )
  {
    START_LOG ("build_sf()", "PetscDMWrapper");

    const DofMap & dof_map = system.get_dof_map();

    const std::vector<dof_id_type> & send_list = dof_map.get_send_list();

    // Number of ghost dofs that send information to this processor
    const PetscInt n_leaves = cast_int<PetscInt>(send_list.size());

    // Number of local dofs, including ghosts dofs
    const PetscInt n_roots = dof_map.n_local_dofs() + n_leaves;

    // This is the vector of dof indices coming from other processors
    // We need to give this to the PetscSF
    // We'll be extra paranoid about this ugly double cast
    static_assert(sizeof(PetscInt) == sizeof(dof_id_type),"PetscInt is not a dof_id_type!");
    PetscInt * local_dofs = reinterpret_cast<PetscInt *>(const_cast<dof_id_type *>(send_list.data()));

    // This is the vector of PetscSFNode's for the local_dofs.
    // For each entry in local_dof, we have to supply the rank from which
    // that dof stems and its local index on that rank.
    // PETSc documentation here:
    // http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/PetscSF/PetscSFNode.html
    std::vector<PetscSFNode> sf_nodes(send_list.size());

    for( unsigned int i = 0; i < send_list.size(); i++ )
      {
        dof_id_type incoming_dof = send_list[i];

        const processor_id_type rank = dof_map.dof_owner(incoming_dof);

        // Dofs are sorted and continuous on the processor so local index
        // is counted up from the first dof on the processor.
        PetscInt index = incoming_dof - dof_map.first_dof(rank);

        sf_nodes[i].rank  = rank; /* Rank of owner */
        sf_nodes[i].index = index;/* Index of dof on rank */
      }

    PetscSFNode * remote_dofs = sf_nodes.data();

    PetscErrorCode ierr;
    ierr = PetscSFCreate(system.comm().get(), &star_forest);LIBMESH_CHKERR(ierr);

    // TODO: We should create pointers to arrays so we don't have to copy
    //       and then can use PETSC_OWN_POINTER where PETSc will take ownership
    //       and delete the memory for us. But then we'd have to use PetscMalloc.
    ierr = PetscSFSetGraph(star_forest,
                           n_roots,
                           n_leaves,
                           local_dofs,
                           PETSC_COPY_VALUES,
                           remote_dofs,
                           PETSC_COPY_VALUES);
    LIBMESH_CHKERR(ierr);

    STOP_LOG ("build_sf()", "PetscDMWrapper");
  }

  void PetscDMWrapper::set_point_range_in_section (const System & system,
                                                   PetscSection & section,
                                                   std::unordered_map<dof_id_type,dof_id_type> & node_map,
                                                   std::unordered_map<dof_id_type,dof_id_type> & elem_map,
                                                   std::map<dof_id_type,unsigned int> & scalar_map)
  {
    // We're expecting this to be empty coming in
    libmesh_assert(node_map.empty());

    // We need to count up the number of active "points" on this processor.
    // Nominally, a "point" in PETSc parlance is a geometric object that can
    // hold DoFs, i.e node, edge, face, elem. Since we handle the mesh and are only
    // interested in solvers, then the only thing PETSc needs is a unique *local* number
    // for each "point" that has active DoFs; note however this local numbering
    // we construct must be continuous.
    //
    // In libMesh, for most finite elements, we just associate those DoFs with the
    // geometric nodes. So can we loop over the nodes on this processor and check
    // if any of the fields are have active DoFs on that node.
    // If so, then we tell PETSc about that "point". At this stage, we just need
    // to count up how many active "points" we have and cache the local number to global id
    // mapping.


    // These will be our local counters. pstart should always be zero.
    // pend will track our local "point" count.
    // If we're on a processor who coarsened the mesh to have no local elements,
    // we should make an empty PetscSection. An empty PetscSection is specified
    // by passing [0,0) to the PetscSectionSetChart call at the end. So, if we
    // have nothing on this processor, these are the correct values to pass to
    // PETSc.
    dof_id_type pstart = 0;
    dof_id_type pend = 0;

    const MeshBase & mesh = system.get_mesh();

    const DofMap & dof_map = system.get_dof_map();

    // If we don't have any local dofs, then there's nothing to tell to the PetscSection
    if (dof_map.n_local_dofs() > 0)
      {
        // Conservative estimate of space needed so we don't thrash
        node_map.reserve(mesh.n_local_nodes());
        elem_map.reserve(mesh.n_active_local_elem());

        // We loop over active elements and then cache the global/local node mapping to make sure
        // we only count active nodes. For example, if we're calling this function and we're
        // not the finest level in the Mesh tree, we don't want to include nodes of child elements
        // that aren't active on this level.
        for (const auto & elem : mesh.active_local_element_ptr_range())
          {
            for (unsigned int n = 0; n < elem->n_nodes(); n++)
              {
                // get the global id number of local node n
                const Node & node = elem->node_ref(n);

                // Only register nodes with the PetscSection if they have dofs that belong to
                // this processor. Even though we're active local elements, the dofs associated
                // with the node may belong to a different processor. The processor who owns
                // those dofs will register that node with the PetscSection on that processor.
                std::vector<dof_id_type> node_dof_indices;
                dof_map.dof_indices( &node, node_dof_indices );
                if( !node_dof_indices.empty() && dof_map.local_index(node_dof_indices[0]) )
                  {
#ifndef NDEBUG
                    // We're assuming that if the first dof for this node belongs to this processor,
                    // then all of them do.
                    for( auto dof : node_dof_indices )
                      libmesh_assert(dof_map.local_index(dof));
#endif
                    // Cache the global/local mapping if we haven't already
                    // Then increment our local count
                    dof_id_type node_id = node.id();
                    if( node_map.count(node_id) == 0 )
                      {
                        node_map.insert(std::make_pair(node_id,pend));
                        pend++;
                      }
                  }
              }

            // Some finite elements, e.g. Hierarchic, associate element interior DoFs with the element
            // rather than the node (since we ought to be able to use Hierachic elements on a QUAD4,
            // which has no interior node). Thus, we also need to check element interiors for DoFs
            // as well and, if the finite element has them, we also need to count the Elem in our
            // "point" accounting.
            if( elem->n_dofs(system.number()) > 0 )
              {
                dof_id_type elem_id = elem->id();
                elem_map.insert(std::make_pair(elem_id,pend));
                pend++;
              }
          }

        // SCALAR dofs live on the "last" processor, so only work there if there are any
        if( dof_map.n_SCALAR_dofs() > 0 && (system.processor_id() == (system.n_processors()-1)) )
          {
            // Loop through all the variables and cache the scalar ones. We cache the
            // SCALAR variable index along with the local point to make it easier when
            // we have to register dofs with the PetscSection
            for( unsigned int v = 0; v < system.n_vars(); v++ )
              {
                if( system.variable(v).type().family == SCALAR )
                  {
                    scalar_map.insert(std::make_pair(pend,v));
                    pend++;
                  }
              }
          }

      }

    PetscErrorCode ierr = PetscSectionSetChart(section, pstart, pend);
    LIBMESH_CHKERR(ierr);
  }

  void PetscDMWrapper::add_dofs_to_section (const System & system,
                                            PetscSection & section,
                                            const std::unordered_map<dof_id_type,dof_id_type> & node_map,
                                            const std::unordered_map<dof_id_type,dof_id_type> & elem_map,
                                            const std::map<dof_id_type,unsigned int> & scalar_map)
  {
    const MeshBase & mesh = system.get_mesh();

    PetscErrorCode ierr;

    // Now we go through and add dof information for each "point".
    //
    // In libMesh, for most finite elements, we just associate those DoFs with the
    // geometric nodes. So can we loop over the nodes we cached in the node_map and
    // the DoFs for each field for that node. We need to give PETSc the local id
    // we built up in the node map.
    for (const auto & nmap : node_map )
      {
        const dof_id_type global_node_id = nmap.first;
        const dof_id_type local_node_id = nmap.second;

        libmesh_assert( mesh.query_node_ptr(global_node_id) );

        const Node & node = mesh.node_ref(global_node_id);

        this->add_dofs_helper(system,node,local_node_id,section);
      }

    // Some finite element types associate dofs with the element. So now we go through
    // any of those with the Elem as the point we add to the PetscSection with accompanying
    // dofs
    for (const auto & emap : elem_map )
      {
        const dof_id_type global_elem_id = emap.first;
        const dof_id_type local_elem_id = emap.second;

        libmesh_assert( mesh.query_elem_ptr(global_elem_id) );

        const Elem & elem = mesh.elem_ref(global_elem_id);

        this->add_dofs_helper(system,elem,local_elem_id,section);
      }

    // Now add any SCALAR dofs to the PetscSection
    // SCALAR dofs live on the "last" processor, so only work there if there are any
    if (system.processor_id() == (system.n_processors()-1))
      {
        for (const auto & smap : scalar_map )
          {
            const dof_id_type local_id = smap.first;
            const unsigned int scalar_var = smap.second;

            // The number of SCALAR dofs comes from the variable order
            const int n_dofs = system.variable(scalar_var).type().order.get_order();

            ierr = PetscSectionSetFieldDof( section, local_id, scalar_var, n_dofs );
            LIBMESH_CHKERR(ierr);

            // In the SCALAR case, there are no other variables associate with the "point"
            // the total number of dofs on the point is the same as that for the field
            ierr = PetscSectionSetDof( section, local_id, n_dofs );
            LIBMESH_CHKERR(ierr);
          }
      }

  }

  void PetscDMWrapper::add_dofs_helper (const System & system,
                                        const DofObject & dof_object,
                                        dof_id_type local_id,
                                        PetscSection & section)
  {
    unsigned int total_n_dofs_at_dofobject = 0;

    // We are assuming variables are also numbered 0 to n_vars()-1
    for( unsigned int v = 0; v < system.n_vars(); v++ )
      {
        unsigned int n_dofs_at_dofobject = dof_object.n_dofs(system.number(), v);

        if( n_dofs_at_dofobject > 0 )
          {
            PetscErrorCode ierr = PetscSectionSetFieldDof( section,
                                                           local_id,
                                                           v,
                                                           n_dofs_at_dofobject );

            LIBMESH_CHKERR(ierr);

            total_n_dofs_at_dofobject += n_dofs_at_dofobject;
          }
      }

    libmesh_assert_equal_to(total_n_dofs_at_dofobject, dof_object.n_dofs(system.number()));

    PetscErrorCode ierr =
      PetscSectionSetDof( section, local_id, total_n_dofs_at_dofobject );
    LIBMESH_CHKERR(ierr);
  }


  dof_id_type PetscDMWrapper::check_section_n_dofs( PetscSection & section )
  {
    PetscInt n_local_dofs = 0;

    // Grap the starting and ending points from the section
    PetscInt pstart, pend;
    PetscErrorCode ierr = PetscSectionGetChart(section, &pstart, &pend);
    LIBMESH_CHKERR(ierr);

    // Count up the n_dofs for each point from the section
    for( PetscInt p = pstart; p < pend; p++ )
      {
        PetscInt n_dofs;
        ierr = PetscSectionGetDof(section,p,&n_dofs);LIBMESH_CHKERR(ierr);
        n_local_dofs += n_dofs;
      }

    static_assert(sizeof(PetscInt) == sizeof(dof_id_type),"PetscInt is not a dof_id_type!");
    return n_local_dofs;
  }

  void PetscDMWrapper::init_dm_data(unsigned int n_levels, const Parallel::Communicator & comm)
  {
    _dms.resize(n_levels);
    _sections.resize(n_levels);
    _star_forests.resize(n_levels);
    _ctx_vec.resize(n_levels);
    _pmtx_vec.resize(n_levels);
    _subpmtx_vec.resize(n_levels);
    _vec_vec.resize(n_levels);
    _mesh_dof_sizes.resize(n_levels);
    _mesh_dof_loc_sizes.resize(n_levels);

    for( unsigned int i = 0; i < n_levels; i++ )
      {
        _dms[i] = libmesh_make_unique<DM>();
        _sections[i] = libmesh_make_unique<PetscSection>();
        _star_forests[i] = libmesh_make_unique<PetscSF>();
        _ctx_vec[i] = libmesh_make_unique<PetscDMContext>();
        _pmtx_vec[i]= libmesh_make_unique<PetscMatrix<Number>>(comm);
        _subpmtx_vec[i]= libmesh_make_unique<PetscMatrix<Number>>(comm);
        _vec_vec[i] = libmesh_make_unique<PetscVector<Number>>(comm);
      }
  }

} // end namespace libMesh

#endif // #if LIBMESH_ENABLE_AMR && LIBMESH_HAVE_METAPHYSICL
#endif // PETSC_VERSION
#endif // LIBMESH_HAVE_PETSC
