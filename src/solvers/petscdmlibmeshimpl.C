// The libMesh Finite Element Library.
// Copyright (C) 2002-2024 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#include "libmesh/petsc_macro.h"

#ifdef LIBMESH_HAVE_PETSC

#include "libmesh/ignore_warnings.h"

// PETSc includes
#if !PETSC_VERSION_LESS_THAN(3,6,0)
# include <petsc/private/dmimpl.h>
#else
# include <petsc-private/dmimpl.h>
#endif

#include "libmesh/restore_warnings.h"

// Local Includes
#include "libmesh/libmesh_common.h"
#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/petsc_nonlinear_solver.h"
#include "libmesh/petsc_linear_solver.h"
#include "libmesh/petsc_vector.h"
#include "libmesh/petsc_matrix.h"
#include "libmesh/petscdmlibmesh.h"
#include "libmesh/dof_map.h"
#include "libmesh/preconditioner.h"
#include "libmesh/elem.h"
#include "libmesh/parallel.h"


using namespace libMesh;


#define DMLIBMESH_NO_DECOMPOSITION     0
#define DMLIBMESH_FIELD_DECOMPOSITION  1
#define DMLIBMESH_DOMAIN_DECOMPOSITION 2

#define DMLIBMESH_NO_EMBEDDING         0
#define DMLIBMESH_FIELD_EMBEDDING      1
#define DMLIBMESH_DOMAIN_EMBEDDING     2

struct DM_libMesh
{
  NonlinearImplicitSystem * sys;
  std::map<std::string, unsigned int> * varids;
  std::map<unsigned int, std::string> * varnames;
  std::map<std::string, unsigned int> * blockids;
  std::map<unsigned int, std::string> * blocknames;
  unsigned int decomposition_type;
  std::vector<std::set<unsigned int>> * decomposition;
  unsigned int embedding_type;
  IS embedding;
  unsigned int vec_count;
};

struct DMVec_libMesh {
  std::string label;
};

#undef  __FUNCT__
#define __FUNCT__ "DMlibMeshGetVec_Private"
PetscErrorCode DMlibMeshGetVec_Private(DM /*dm*/, const char * /*name*/, Vec *)
{
  PetscFunctionBegin;

  PetscFunctionReturn(LIBMESH_PETSC_SUCCESS);
}



PetscErrorCode DMlibMeshSetUpName_Private(DM dm);

#undef  __FUNCT__
#define __FUNCT__ "DMlibMeshSetSystem_libMesh"
PetscErrorCode DMlibMeshSetSystem_libMesh(DM dm, NonlinearImplicitSystem & sys)
{
  const Parallel::Communicator & comm = sys.comm();

  PetscErrorCode ierr;
  PetscFunctionBegin;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  PetscBool islibmesh;
  ierr = PetscObjectTypeCompare((PetscObject)dm, DMLIBMESH,&islibmesh);
  CHKERRQ(ierr);
  if (!islibmesh) LIBMESH_SETERRQ2(((PetscObject)dm)->comm, PETSC_ERR_ARG_WRONG, "Got DM of type %s, not of type %s", ((PetscObject)dm)->type_name, DMLIBMESH);

  if (dm->setupcalled) SETERRQ(((PetscObject)dm)->comm, PETSC_ERR_ARG_WRONGSTATE, "Cannot reset the libMesh system after DM has been set up.");
  DM_libMesh * dlm = (DM_libMesh *)(dm->data);
  dlm->sys =&sys;
  /* Initially populate the sets of active blockids and varids using all of the
     existing blocks/variables (only variables are supported at the moment). */
  DofMap & dofmap = dlm->sys->get_dof_map();
  dlm->varids->clear();
  dlm->varnames->clear();
  for (auto v : make_range(dofmap.n_variables())) {
    std::string vname = dofmap.variable(v).name();
    dlm->varids->insert(std::pair<std::string,unsigned int>(vname,v));
    dlm->varnames->insert(std::pair<unsigned int,std::string>(v,vname));
  }
  const MeshBase & mesh = dlm->sys->get_mesh();
  dlm->blockids->clear();
  dlm->blocknames->clear();
  std::set<subdomain_id_type> blocks;
  /* The following effectively is a verbatim copy of MeshBase::n_subdomains(). */
  // This requires an inspection on every processor
  libmesh_parallel_only(mesh.comm());
  for (const auto & elem : mesh.active_element_ptr_range())
    blocks.insert(elem->subdomain_id());
  // Some subdomains may only live on other processors
  comm.set_union(blocks);

  std::set<subdomain_id_type>::iterator bit = blocks.begin();
  std::set<subdomain_id_type>::iterator bend = blocks.end();
  if (bit == bend) SETERRQ(((PetscObject)dm)->comm, PETSC_ERR_PLIB, "No mesh blocks found.");

  for (; bit != bend; ++bit) {
    subdomain_id_type bid = *bit;
    std::string bname = mesh.subdomain_name(bid);
    if (!bname.length()) {
      /* Block names are currently implemented for Exodus II meshes
         only, so we might have to make up our own block names and
         maintain our own mapping of block ids to names.
      */
      std::ostringstream ss;
      ss << "dm" << bid;
      bname = ss.str();
    }
    dlm->blockids->insert(std::pair<std::string,unsigned int>(bname,bid));
    dlm->blocknames->insert(std::pair<unsigned int,std::string>(bid,bname));
  }
  ierr = DMlibMeshSetUpName_Private(dm); CHKERRQ(ierr);
  PetscFunctionReturn(LIBMESH_PETSC_SUCCESS);
}

#undef  __FUNCT__
#define __FUNCT__ "DMlibMeshGetSystem_libMesh"
PetscErrorCode DMlibMeshGetSystem_libMesh(DM dm, NonlinearImplicitSystem *& sys)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  PetscBool islibmesh;
  ierr = PetscObjectTypeCompare((PetscObject)dm, DMLIBMESH,&islibmesh);CHKERRQ(ierr);
  if (!islibmesh) LIBMESH_SETERRQ2(((PetscObject)dm)->comm, PETSC_ERR_ARG_WRONG, "Got DM of type %s, not of type %s", ((PetscObject)dm)->type_name, DMLIBMESH);
  DM_libMesh * dlm = (DM_libMesh *)(dm->data);
  sys = dlm->sys;
  PetscFunctionReturn(LIBMESH_PETSC_SUCCESS);
}


#undef  __FUNCT__
#define __FUNCT__ "DMlibMeshGetBlocks"
PetscErrorCode DMlibMeshGetBlocks(DM dm, PetscInt * n, char *** blocknames)
{
  PetscErrorCode ierr;
  PetscInt i;
  PetscFunctionBegin;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  PetscBool islibmesh;
  ierr = PetscObjectTypeCompare((PetscObject)dm, DMLIBMESH,&islibmesh);
  CHKERRQ(ierr);
  if (!islibmesh) LIBMESH_SETERRQ2(((PetscObject)dm)->comm, PETSC_ERR_ARG_WRONG, "Got DM of type %s, not of type %s", ((PetscObject)dm)->type_name, DMLIBMESH);
  DM_libMesh * dlm = (DM_libMesh *)(dm->data);
#if PETSC_RELEASE_GREATER_EQUALS(3,20,0)
  PetscAssertPointer(n,2);
#else
  PetscValidPointer(n,2);
#endif
  *n = cast_int<unsigned int>(dlm->blockids->size());
  if (!blocknames) PetscFunctionReturn(LIBMESH_PETSC_SUCCESS);
  ierr = PetscMalloc(*n*sizeof(char *), blocknames); CHKERRQ(ierr);
  i = 0;
  for (const auto & pr : *(dlm->blockids))
    {
      ierr = PetscStrallocpy(pr.first.c_str(), *blocknames+i); CHKERRQ(ierr);
      ++i;
    }
  PetscFunctionReturn(LIBMESH_PETSC_SUCCESS);
}

#undef  __FUNCT__
#define __FUNCT__ "DMlibMeshGetVariables"
PetscErrorCode DMlibMeshGetVariables(DM dm, PetscInt * n, char *** varnames)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  PetscBool islibmesh;
  PetscInt i;
  ierr = PetscObjectTypeCompare((PetscObject)dm, DMLIBMESH,&islibmesh);
  CHKERRQ(ierr);
  if (!islibmesh) LIBMESH_SETERRQ2(((PetscObject)dm)->comm, PETSC_ERR_ARG_WRONG, "Got DM of type %s, not of type %s", ((PetscObject)dm)->type_name, DMLIBMESH);
  DM_libMesh * dlm = (DM_libMesh *)(dm->data);
#if PETSC_RELEASE_GREATER_EQUALS(3,20,0)
  PetscAssertPointer(n,2);
#else
  PetscValidPointer(n,2);
#endif
  *n = cast_int<unsigned int>(dlm->varids->size());
  if (!varnames) PetscFunctionReturn(LIBMESH_PETSC_SUCCESS);
  ierr = PetscMalloc(*n*sizeof(char *), varnames); CHKERRQ(ierr);
  i = 0;
  for (const auto & pr : *(dlm->varids))
    {
      ierr = PetscStrallocpy(pr.first.c_str(), *varnames+i); CHKERRQ(ierr);
      ++i;
    }
  PetscFunctionReturn(LIBMESH_PETSC_SUCCESS);
}

#undef  __FUNCT__
#define __FUNCT__ "DMlibMeshSetUpName_Private"
PetscErrorCode DMlibMeshSetUpName_Private(DM dm)
{
  DM_libMesh * dlm = (DM_libMesh *)dm->data;
  PetscErrorCode ierr;
  PetscFunctionBegin;
  std::string name = dlm->sys->name();
  std::map<unsigned int, std::string> * dnames = LIBMESH_PETSC_NULLPTR,
                                      * enames = LIBMESH_PETSC_NULLPTR;
  if (dlm->decomposition_type == DMLIBMESH_FIELD_DECOMPOSITION) {
    name += ":dec:var:";
    dnames = dlm->varnames;
  }
  if (dlm->decomposition_type == DMLIBMESH_DOMAIN_DECOMPOSITION) {
    name += ":dec:block:";
    dnames = dlm->blocknames;
  }
  if (dnames) {
    for (auto decomp : *dlm->decomposition) {
      for (std::set<unsigned int>::iterator dit_begin = decomp.begin(),
                                            dit = dit_begin,
                                            dit_end   = decomp.end();
           dit != dit_end; ++dit) {
        unsigned int id = *dit;
        if (dit != dit_begin)
          name += ",";
        name += (*dnames)[id];
      }
      name += ";";
    }
  }
  if (dlm->embedding_type == DMLIBMESH_FIELD_EMBEDDING) {
    name += ":emb:var:";
    enames = dlm->varnames;
  }
  if (dlm->embedding_type == DMLIBMESH_DOMAIN_EMBEDDING) {
    name += ":emb:block:";
    enames = dlm->blocknames;
  }
  if (enames) {
    for (auto eit     = enames->begin(),
              eit_end = enames->end(); eit != eit_end; ++eit) {
      std::string & ename = eit->second;
      if (eit != enames->begin())
        name += ",";
      name += ename;
    }
    name += ";";
  }
  ierr = PetscObjectSetName((PetscObject)dm, name.c_str()); CHKERRQ(ierr);
  PetscFunctionReturn(LIBMESH_PETSC_SUCCESS);
}


#undef __FUNCT__
#define __FUNCT__ "DMCreateFieldDecomposition_libMesh"
static PetscErrorCode  DMCreateFieldDecomposition_libMesh(DM dm, PetscInt * len, char *** namelist, IS ** islist, DM ** dmlist)
{
  PetscFunctionBegin;
  PetscErrorCode ierr;
  DM_libMesh     * dlm = (DM_libMesh *)(dm->data);
  NonlinearImplicitSystem * sys = dlm->sys;
  IS emb;
  if (dlm->decomposition_type != DMLIBMESH_FIELD_DECOMPOSITION) PetscFunctionReturn(LIBMESH_PETSC_SUCCESS);

  *len = cast_int<unsigned int>(dlm->decomposition->size());
  if (namelist) {ierr = PetscMalloc(*len*sizeof(char *), namelist);  CHKERRQ(ierr);}
  if (islist)   {ierr = PetscMalloc(*len*sizeof(IS),    islist);    CHKERRQ(ierr);}
  if (dmlist)   {ierr = PetscMalloc(*len*sizeof(DM),    dmlist);    CHKERRQ(ierr);}
  DofMap & dofmap = dlm->sys->get_dof_map();
  for (auto d : index_range(*dlm->decomposition)) {
    std::set<numeric_index_type>         dindices;
    std::string                          dname;
    std::map<std::string, unsigned int>  dvarids;
    std::map<unsigned int, std::string>  dvarnames;
    unsigned int                         dvcount = 0;
    for (const auto & v : (*dlm->decomposition)[d]) {
      std::string vname = (*dlm->varnames)[v];
      dvarids.insert(std::pair<std::string, unsigned int>(vname,v));
      dvarnames.insert(std::pair<unsigned int,std::string>(v,vname));
      if (!dvcount) dname = vname;
      else   dname += "_" + vname;
      ++dvcount;
      if (!islist) continue;
      // Iterate only over this DM's blocks.
      for (const auto & pr : *(dlm->blockids)) {
        const subdomain_id_type sbd_id = cast_int<subdomain_id_type>(pr.second);
        for (const auto & elem :
               as_range(sys->get_mesh().active_local_subdomain_elements_begin(sbd_id),
                        sys->get_mesh().active_local_subdomain_elements_end(sbd_id))) {
          //unsigned int e_subdomain = elem->subdomain_id();
          std::vector<numeric_index_type> evindices;
          // Get the degree of freedom indices for the given variable off the current element.
          dofmap.dof_indices(elem, evindices, v);
          for (numeric_index_type dof : evindices) {
            if (dof >= dofmap.first_dof() && dof < dofmap.end_dof()) // might want to use variable_first/last_local_dof instead
              dindices.insert(dof);
          }
        }
      }
    }
    if (namelist) {
      ierr = PetscStrallocpy(dname.c_str(),(*namelist)+d);            CHKERRQ(ierr);
    }
    if (islist) {
      IS dis;
      PetscInt * darray;
      ierr = PetscMalloc(sizeof(PetscInt)*dindices.size(), &darray); CHKERRQ(ierr);
      numeric_index_type i = 0;
      for (const auto & id : dindices) {
        darray[i] = id;
        ++i;
      }
      ierr = ISCreateGeneral(((PetscObject)dm)->comm,
                             cast_int<PetscInt>(dindices.size()),
                             darray, PETSC_OWN_POINTER, &dis);
      CHKERRQ(ierr);
      if (dlm->embedding) {
        /* Create a relative embedding into the parent's index space. */
        ierr = ISEmbed(dis,dlm->embedding, PETSC_TRUE, &emb); CHKERRQ(ierr);
        PetscInt elen, dlen;
        ierr = ISGetLocalSize(emb, &elen); CHKERRQ(ierr);
        ierr = ISGetLocalSize(dis, &dlen); CHKERRQ(ierr);
        if (elen != dlen) LIBMESH_SETERRQ1(((PetscObject)dm)->comm, PETSC_ERR_PLIB, "Failed to embed subdomain %zu", d);
        ierr = ISDestroy(&dis); CHKERRQ(ierr);
        dis = emb;
      }
      else {
        emb = dis;
      }
      (*islist)[d] = dis;
    }
    if (dmlist) {
      DM ddm;
      ierr = DMCreate(((PetscObject)dm)->comm, &ddm); CHKERRQ(ierr);
      ierr = DMSetType(ddm, DMLIBMESH);               CHKERRQ(ierr);
      DM_libMesh * ddlm = (DM_libMesh *)(ddm->data);
      ddlm->sys = dlm->sys;
      /* copy over the block ids and names */
      *ddlm->blockids = *dlm->blockids;
      *ddlm->blocknames = *dlm->blocknames;
      /* set the vars from the d-th part of the decomposition. */
      *ddlm->varids     = dvarids;
      *ddlm->varnames   = dvarnames;
      ierr = PetscObjectReference((PetscObject)emb); CHKERRQ(ierr);
      ddlm->embedding = emb;
      ddlm->embedding_type = DMLIBMESH_FIELD_EMBEDDING;

      ierr = DMlibMeshSetUpName_Private(ddm); CHKERRQ(ierr);
      ierr = DMSetFromOptions(ddm);           CHKERRQ(ierr);
      (*dmlist)[d] = ddm;
    }
  }
  PetscFunctionReturn(LIBMESH_PETSC_SUCCESS);
}

#undef __FUNCT__
#define __FUNCT__ "DMCreateDomainDecomposition_libMesh"
static PetscErrorCode  DMCreateDomainDecomposition_libMesh(DM dm, PetscInt * len, char *** namelist, IS ** innerislist, IS ** outerislist, DM ** dmlist)
{
  PetscFunctionBegin;
  PetscErrorCode ierr;
  DM_libMesh     * dlm = (DM_libMesh *)(dm->data);
  NonlinearImplicitSystem * sys = dlm->sys;
  IS emb;
  if (dlm->decomposition_type != DMLIBMESH_DOMAIN_DECOMPOSITION) PetscFunctionReturn(LIBMESH_PETSC_SUCCESS);
  *len = cast_int<unsigned int>(dlm->decomposition->size());
  if (namelist)      {ierr = PetscMalloc(*len*sizeof(char *), namelist);  CHKERRQ(ierr);}
  if (innerislist)   {ierr = PetscMalloc(*len*sizeof(IS),    innerislist);    CHKERRQ(ierr);}
  if (outerislist)   *outerislist = LIBMESH_PETSC_NULLPTR; /* FIX: allow mesh-based overlap. */
  if (dmlist)        {ierr = PetscMalloc(*len*sizeof(DM),    dmlist);    CHKERRQ(ierr);}
  for (auto d : index_range(*dlm->decomposition)) {
    std::set<numeric_index_type>               dindices;
    std::string                          dname;
    std::map<std::string, unsigned int>  dblockids;
    std::map<unsigned int,std::string>   dblocknames;
    unsigned int                         dbcount = 0;
    for (const auto & b : (*dlm->decomposition)[d]) {
      std::string bname = (*dlm->blocknames)[b];
      dblockids.insert(std::pair<std::string, unsigned int>(bname,b));
      dblocknames.insert(std::pair<unsigned int,std::string>(b,bname));
      if (!dbcount) dname = bname;
      else   dname += "_" + bname;
      ++dbcount;
      if (!innerislist) continue;
      const subdomain_id_type b_sbd_id = cast_int<subdomain_id_type>(b);
      MeshBase::const_element_iterator       el     = sys->get_mesh().active_local_subdomain_elements_begin(b_sbd_id);
      const MeshBase::const_element_iterator end_el = sys->get_mesh().active_local_subdomain_elements_end(b_sbd_id);
      for ( ; el != end_el; ++el) {
        const Elem * elem = *el;
        std::vector<numeric_index_type> evindices;
        // Iterate only over this DM's variables.
        for (const auto & pr : *(dlm->varids)) {
          // Get the degree of freedom indices for the given variable off the current element.
          sys->get_dof_map().dof_indices(elem, evindices, pr.second);
          for (const auto & dof : evindices) {
            if (dof >= sys->get_dof_map().first_dof() && dof < sys->get_dof_map().end_dof()) // might want to use variable_first/last_local_dof instead
              dindices.insert(dof);
          }
        }
      }
    }
    if (namelist) {
      ierr = PetscStrallocpy(dname.c_str(),(*namelist)+d);            CHKERRQ(ierr);
    }
    if (innerislist) {
      PetscInt * darray;
      IS dis;
      ierr = PetscMalloc(sizeof(PetscInt)*dindices.size(), &darray); CHKERRQ(ierr);
      numeric_index_type i = 0;
      for (const auto & id : dindices) {
        darray[i] = id;
        ++i;
      }
      ierr = ISCreateGeneral(((PetscObject)dm)->comm,
                             cast_int<PetscInt>(dindices.size()),
                             darray, PETSC_OWN_POINTER, &dis);
      CHKERRQ(ierr);
      if (dlm->embedding) {
        /* Create a relative embedding into the parent's index space. */
        ierr = ISEmbed(dis,dlm->embedding, PETSC_TRUE, &emb); CHKERRQ(ierr);
        PetscInt elen, dlen;
        ierr = ISGetLocalSize(emb, &elen); CHKERRQ(ierr);
        ierr = ISGetLocalSize(dis, &dlen);  CHKERRQ(ierr);
        if (elen != dlen) LIBMESH_SETERRQ1(((PetscObject)dm)->comm, PETSC_ERR_PLIB, "Failed to embed field %zu" , d);
        ierr = ISDestroy(&dis); CHKERRQ(ierr);
        dis = emb;
      }
      else {
        emb = dis;
      }
      if (innerislist) {
        ierr = PetscObjectReference((PetscObject)dis); CHKERRQ(ierr);
        (*innerislist)[d] = dis;
      }
      ierr = ISDestroy(&dis); CHKERRQ(ierr);
    }
    if (dmlist) {
      DM ddm;
      ierr = DMCreate(((PetscObject)dm)->comm, &ddm); CHKERRQ(ierr);
      ierr = DMSetType(ddm, DMLIBMESH);               CHKERRQ(ierr);
      DM_libMesh * ddlm = (DM_libMesh *)(ddm->data);
      ddlm->sys = dlm->sys;
      /* copy over the varids and varnames */
      *ddlm->varids    = *dlm->varids;
      *ddlm->varnames  = *dlm->varnames;
      /* set the blocks from the d-th part of the decomposition. */
      *ddlm->blockids    = dblockids;
      *ddlm->blocknames  = dblocknames;
      ierr = PetscObjectReference((PetscObject)emb); CHKERRQ(ierr);
      ddlm->embedding = emb;
      ddlm->embedding_type = DMLIBMESH_DOMAIN_EMBEDDING;

      ierr = DMlibMeshSetUpName_Private(ddm); CHKERRQ(ierr);
      ierr = DMSetFromOptions(ddm);           CHKERRQ(ierr);
      (*dmlist)[d] = ddm;
    }
  }
  PetscFunctionReturn(LIBMESH_PETSC_SUCCESS);
}


#undef  __FUNCT__
#define __FUNCT__ "DMlibMeshCreateFieldDecompositionDM"
PetscErrorCode DMlibMeshCreateFieldDecompositionDM(DM dm, PetscInt dnumber, PetscInt * dsizes, char *** dvarlists, DM * ddm)
{
  PetscErrorCode ierr;
  PetscBool islibmesh;
  PetscFunctionBegin;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  ierr = PetscObjectTypeCompare((PetscObject)dm, DMLIBMESH,&islibmesh);
  CHKERRQ(ierr);
  if (!islibmesh) LIBMESH_SETERRQ2(((PetscObject)dm)->comm, PETSC_ERR_ARG_WRONG, "Got DM of type %s, not of type %s", ((PetscObject)dm)->type_name, DMLIBMESH);
  if (dnumber < 0) LIBMESH_SETERRQ1(((PetscObject)dm)->comm, PETSC_ERR_ARG_WRONG, "Negative number %" LIBMESH_PETSCINT_FMT " of decomposition parts", dnumber);
#if PETSC_RELEASE_GREATER_EQUALS(3,20,0)
  PetscAssertPointer(ddm,5);
#else
  PetscValidPointer(ddm,5);
#endif
  DM_libMesh * dlm = (DM_libMesh *)(dm->data);
  ierr = DMCreate(((PetscObject)dm)->comm, ddm); CHKERRQ(ierr);
  ierr = DMSetType(*ddm, DMLIBMESH);             CHKERRQ(ierr);
  DM_libMesh * ddlm = (DM_libMesh *)((*ddm)->data);
  ddlm->sys = dlm->sys;
  ddlm->varids = dlm->varids;
  ddlm->varnames = dlm->varnames;
  ddlm->blockids = dlm->blockids;
  ddlm->blocknames = dlm->blocknames;
  ddlm->decomposition = new(std::vector<std::set<unsigned int>>);
  ddlm->decomposition_type = DMLIBMESH_FIELD_DECOMPOSITION;
  if (dnumber) {
    for (PetscInt d = 0; d < dnumber; ++d) {
      if (dsizes[d] < 0) LIBMESH_SETERRQ2(((PetscObject)dm)->comm, PETSC_ERR_ARG_WRONG, "Negative size %" LIBMESH_PETSCINT_FMT " of decomposition part %" LIBMESH_PETSCINT_FMT, dsizes[d],d);
      ddlm->decomposition->push_back(std::set<unsigned int>());
      for (PetscInt v = 0; v < dsizes[d]; ++v) {
        std::string vname(dvarlists[d][v]);
        std::map<std::string, unsigned int>::const_iterator vit = dlm->varids->find(vname);
        if (vit == dlm->varids->end())
          LIBMESH_SETERRQ3(((PetscObject)dm)->comm, PETSC_ERR_ARG_WRONG, "Variable %" LIBMESH_PETSCINT_FMT " on the %" LIBMESH_PETSCINT_FMT "-th list with name %s is not owned by this DM", v, d, dvarlists[d][v]);
        unsigned int vid = vit->second;
        (*ddlm->decomposition)[d].insert(vid);
      }
    }
  }
  else { // Empty splits indicate default: split all variables with one per split.
    PetscInt d = 0;
    for (const auto & pr : (*ddlm->varids)) {
      ddlm->decomposition->push_back(std::set<unsigned int>());
      unsigned int vid = pr.second;
      (*ddlm->decomposition)[d].insert(vid);
      ++d;
    }
  }
  ierr = DMlibMeshSetUpName_Private(*ddm); CHKERRQ(ierr);
  ierr = DMSetFromOptions(*ddm);           CHKERRQ(ierr);
  ierr = DMSetUp(*ddm);                    CHKERRQ(ierr);
  PetscFunctionReturn(LIBMESH_PETSC_SUCCESS);
}

#undef  __FUNCT__
#define __FUNCT__ "DMlibMeshCreateDomainDecompositionDM"
PetscErrorCode DMlibMeshCreateDomainDecompositionDM(DM dm, PetscInt dnumber, PetscInt * dsizes, char *** dblocklists, DM * ddm)
{
  PetscErrorCode ierr;
  PetscBool islibmesh;
  PetscFunctionBegin;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  ierr = PetscObjectTypeCompare((PetscObject)dm, DMLIBMESH,&islibmesh);
  CHKERRQ(ierr);
  if (!islibmesh) LIBMESH_SETERRQ2(((PetscObject)dm)->comm, PETSC_ERR_ARG_WRONG, "Got DM of type %s, not of type %s", ((PetscObject)dm)->type_name, DMLIBMESH);
  if (dnumber < 0) LIBMESH_SETERRQ1(((PetscObject)dm)->comm, PETSC_ERR_ARG_WRONG, "Negative number %" LIBMESH_PETSCINT_FMT " of decomposition parts", dnumber);
#if PETSC_RELEASE_GREATER_EQUALS(3,20,0)
  PetscAssertPointer(ddm,5);
#else
  PetscValidPointer(ddm,5);
#endif
  DM_libMesh * dlm = (DM_libMesh *)(dm->data);
  ierr = DMCreate(((PetscObject)dm)->comm, ddm); CHKERRQ(ierr);
  ierr = DMSetType(*ddm, DMLIBMESH);             CHKERRQ(ierr);
  DM_libMesh * ddlm = (DM_libMesh *)((*ddm)->data);
  ddlm->sys = dlm->sys;
  ddlm->varids   = dlm->varids;
  ddlm->varnames = dlm->varnames;
  ddlm->blockids   = dlm->blockids;
  ddlm->blocknames = dlm->blocknames;
  ddlm->decomposition = new(std::vector<std::set<unsigned int>>);
  ddlm->decomposition_type = DMLIBMESH_DOMAIN_DECOMPOSITION;
  if (dnumber) {
    for (PetscInt d = 0; d < dnumber; ++d) {
      if (dsizes[d] < 0) LIBMESH_SETERRQ2(((PetscObject)dm)->comm, PETSC_ERR_ARG_WRONG, "Negative size %" LIBMESH_PETSCINT_FMT " of decomposition part %" LIBMESH_PETSCINT_FMT, dsizes[d],d);
      ddlm->decomposition->push_back(std::set<unsigned int>());
      for (PetscInt b = 0; b < dsizes[d]; ++b) {
        std::string bname(dblocklists[d][b]);
        std::map<std::string, unsigned int>::const_iterator bit = dlm->blockids->find(bname);
        if (bit == dlm->blockids->end())
          LIBMESH_SETERRQ3(((PetscObject)dm)->comm, PETSC_ERR_ARG_WRONG, "Block %" LIBMESH_PETSCINT_FMT " on the %" LIBMESH_PETSCINT_FMT "-th list with name %s is not owned by this DM", b, d, dblocklists[d][b]);
        unsigned int bid = bit->second;
        (*ddlm->decomposition)[d].insert(bid);
      }
    }
  }
  else { // Empty splits indicate default: split all blocks with one per split.
    PetscInt d = 0;
    for (const auto & pr : (*ddlm->blockids)) {
      ddlm->decomposition->push_back(std::set<unsigned int>());
      unsigned int bid = pr.second;
      (*ddlm->decomposition)[d].insert(bid);
      ++d;
    }
  }
  ierr = DMlibMeshSetUpName_Private(*ddm); CHKERRQ(ierr);
  ierr = DMSetFromOptions(*ddm);           CHKERRQ(ierr);
  ierr = DMSetUp(*ddm);                    CHKERRQ(ierr);
  PetscFunctionReturn(LIBMESH_PETSC_SUCCESS);
}

struct token {
  const  char * s;
  struct token * next;
};



#undef __FUNCT__
#define __FUNCT__ "DMlibMeshFunction"
static PetscErrorCode DMlibMeshFunction(DM dm, Vec x, Vec r)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  libmesh_assert(x);
  libmesh_assert(r);

  NonlinearImplicitSystem * _sys;
  ierr = DMlibMeshGetSystem(dm, _sys);CHKERRQ(ierr);
  NonlinearImplicitSystem & sys = *_sys;
  PetscVector<Number> & X_sys = *cast_ptr<PetscVector<Number> *>(sys.solution.get());
  PetscVector<Number> & R_sys = *cast_ptr<PetscVector<Number> *>(sys.rhs);
  PetscVector<Number> X_global(x, _sys->comm()), R(r, _sys->comm());

  // Use the systems update() to get a good local version of the parallel solution
  X_global.swap(X_sys);
  R.swap(R_sys);

  _sys->get_dof_map().enforce_constraints_exactly(*_sys);
  _sys->update();

  // Swap back
  X_global.swap(X_sys);
  R.swap(R_sys);
  R.zero();

  // if the user has provided both function pointers and objects only the pointer
  // will be used, so catch that as an error
  libmesh_error_msg_if(_sys->nonlinear_solver->residual && _sys->nonlinear_solver->residual_object,
                       "ERROR: cannot specify both a function and object to compute the Residual!");

  libmesh_error_msg_if(_sys->nonlinear_solver->matvec && _sys->nonlinear_solver->residual_and_jacobian_object,
                       "ERROR: cannot specify both a function and object to compute the combined Residual & Jacobian!");

  if (_sys->nonlinear_solver->residual != nullptr)
    _sys->nonlinear_solver->residual(*(_sys->current_local_solution.get()), R, *_sys);

  else if (_sys->nonlinear_solver->residual_object != nullptr)
    _sys->nonlinear_solver->residual_object->residual(*(_sys->current_local_solution.get()), R, *_sys);

  else if (_sys->nonlinear_solver->matvec   != nullptr)
    _sys->nonlinear_solver->matvec(*(_sys->current_local_solution.get()), &R, nullptr, *_sys);

  else if (_sys->nonlinear_solver->residual_and_jacobian_object != nullptr)
    _sys->nonlinear_solver->residual_and_jacobian_object->residual_and_jacobian(*(_sys->current_local_solution.get()), &R, nullptr, *_sys);

  else
    libmesh_error_msg("Error! Unable to compute residual and/or Jacobian!");

  R.close();
  X_global.close();
  PetscFunctionReturn(LIBMESH_PETSC_SUCCESS);
}

#undef __FUNCT__
#define __FUNCT__ "SNESFunction_DMlibMesh"
static PetscErrorCode SNESFunction_DMlibMesh(SNES, Vec x, Vec r, void * ctx)
{
  DM dm = (DM)ctx;
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = DMlibMeshFunction(dm,x,r);CHKERRQ(ierr);
  PetscFunctionReturn(LIBMESH_PETSC_SUCCESS);
}


#undef __FUNCT__
#define __FUNCT__ "DMlibMeshJacobian"
static PetscErrorCode DMlibMeshJacobian(DM dm, Vec x, Mat jac, Mat pc)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  NonlinearImplicitSystem * _sys;
  ierr = DMlibMeshGetSystem(dm, _sys); CHKERRQ(ierr);
  NonlinearImplicitSystem & sys = *_sys;

  PetscMatrix<Number> the_pc(pc,sys.comm());
  PetscMatrix<Number> Jac(jac,sys.comm());
  PetscVector<Number> & X_sys = *cast_ptr<PetscVector<Number> *>(sys.solution.get());
  PetscMatrix<Number> & Jac_sys = *cast_ptr<PetscMatrix<Number> *>(sys.matrix);
  PetscVector<Number> X_global(x, sys.comm());

  // Set the dof maps
  the_pc.attach_dof_map(sys.get_dof_map());
  Jac.attach_dof_map(sys.get_dof_map());

  // Use the systems update() to get a good local version of the parallel solution
  X_global.swap(X_sys);
  Jac.swap(Jac_sys);

  sys.get_dof_map().enforce_constraints_exactly(sys);
  sys.update();

  X_global.swap(X_sys);
  Jac.swap(Jac_sys);

  the_pc.zero();

  // if the user has provided both function pointers and objects only the pointer
  // will be used, so catch that as an error
  libmesh_error_msg_if(sys.nonlinear_solver->jacobian && sys.nonlinear_solver->jacobian_object,
                       "ERROR: cannot specify both a function and object to compute the Jacobian!");

  libmesh_error_msg_if(sys.nonlinear_solver->matvec && sys.nonlinear_solver->residual_and_jacobian_object,
                       "ERROR: cannot specify both a function and object to compute the combined Residual & Jacobian!");

  if (sys.nonlinear_solver->jacobian != nullptr)
    sys.nonlinear_solver->jacobian(*(sys.current_local_solution.get()), the_pc, sys);

  else if (sys.nonlinear_solver->jacobian_object != nullptr)
    sys.nonlinear_solver->jacobian_object->jacobian(*(sys.current_local_solution.get()), the_pc, sys);

  else if (sys.nonlinear_solver->matvec != nullptr)
    sys.nonlinear_solver->matvec(*(sys.current_local_solution.get()), nullptr, &the_pc, sys);

  else if (sys.nonlinear_solver->residual_and_jacobian_object != nullptr)
    sys.nonlinear_solver->residual_and_jacobian_object->residual_and_jacobian(*(sys.current_local_solution.get()), nullptr, &the_pc, sys);

  else
    libmesh_error_msg("Error! Unable to compute residual and/or Jacobian!");

  the_pc.close();
  Jac.close();
  X_global.close();
  PetscFunctionReturn(LIBMESH_PETSC_SUCCESS);
}

#undef  __FUNCT__
#define __FUNCT__ "SNESJacobian_DMlibMesh"
static PetscErrorCode SNESJacobian_DMlibMesh(SNES, Vec x, Mat jac, Mat pc, void * ctx)
{
  DM dm = (DM)ctx;
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = DMlibMeshJacobian(dm,x,jac,pc); CHKERRQ(ierr);
  PetscFunctionReturn(LIBMESH_PETSC_SUCCESS);
}

#undef __FUNCT__
#define __FUNCT__ "DMVariableBounds_libMesh"
static PetscErrorCode DMVariableBounds_libMesh(DM dm, Vec xl, Vec xu)
{
  PetscErrorCode ierr;
  NonlinearImplicitSystem * _sys;
  ierr = DMlibMeshGetSystem(dm, _sys); CHKERRQ(ierr);
  NonlinearImplicitSystem & sys = *_sys;
  PetscVector<Number> XL(xl, sys.comm());
  PetscVector<Number> XU(xu, sys.comm());
  PetscFunctionBegin;
  // Workaround for nonstandard Q suffix warning with quad precision
  const PetscReal petsc_inf = std::numeric_limits<PetscReal>::max() / 4;
  ierr = VecSet(xl, -petsc_inf);CHKERRQ(ierr);
  ierr = VecSet(xu, petsc_inf);CHKERRQ(ierr);
  if (sys.nonlinear_solver->bounds != nullptr)
    sys.nonlinear_solver->bounds(XL,XU,sys);
  else if (sys.nonlinear_solver->bounds_object != nullptr)
    sys.nonlinear_solver->bounds_object->bounds(XL,XU, sys);
  else
    SETERRQ(((PetscObject)dm)->comm, PETSC_ERR_ARG_WRONG, "No bounds calculation in this libMesh object");

  PetscFunctionReturn(LIBMESH_PETSC_SUCCESS);
}


#undef __FUNCT__
#define __FUNCT__ "DMCreateGlobalVector_libMesh"
static PetscErrorCode DMCreateGlobalVector_libMesh(DM dm, Vec *x)
{
  PetscFunctionBegin;
  PetscErrorCode ierr;
  DM_libMesh     * dlm = (DM_libMesh *)(dm->data);
  PetscBool eq;

  ierr = PetscObjectTypeCompare((PetscObject)dm, DMLIBMESH, &eq); CHKERRQ(ierr);

  if (!eq)
    LIBMESH_SETERRQ2(((PetscObject)dm)->comm, PETSC_ERR_ARG_WRONG, "DM of type %s, not of type %s", ((PetscObject)dm)->type_name, DMLIBMESH);

  if (!dlm->sys)
    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONGSTATE, "No libMesh system set for DM_libMesh");

  NumericVector<Number> * nv = (dlm->sys->solution).get();
  PetscVector<Number> *   pv = dynamic_cast<PetscVector<Number> *>(nv);
  Vec                    v  = pv->vec();
  /* Unfortunately, currently this does not produce a ghosted vector, so nonlinear subproblem solves aren't going to be easily available.
     Should work fine for getting vectors out for linear subproblem solvers. */
  if (dlm->embedding) {
    PetscInt n;
    ierr = VecCreate(((PetscObject)v)->comm, x);       CHKERRQ(ierr);
    ierr = ISGetLocalSize(dlm->embedding, &n);         CHKERRQ(ierr);
    ierr = VecSetSizes(*x,n,PETSC_DETERMINE);           CHKERRQ(ierr);
    ierr = VecSetType(*x,((PetscObject)v)->type_name); CHKERRQ(ierr);
    ierr = VecSetFromOptions(*x);                      CHKERRQ(ierr);
    ierr = VecSetUp(*x);                               CHKERRQ(ierr);
  }
  else {
    ierr = VecDuplicate(v,x); CHKERRQ(ierr);
  }

#if PETSC_VERSION_LESS_THAN(3,13,0)
  ierr = PetscObjectCompose((PetscObject)*x,"DM",(PetscObject)dm); CHKERRQ(ierr);
#else
  ierr = VecSetDM(*x, dm);CHKERRQ(ierr);
#endif

  PetscFunctionReturn(LIBMESH_PETSC_SUCCESS);
}




#undef __FUNCT__
#define __FUNCT__ "DMCreateMatrix_libMesh"
static PetscErrorCode DMCreateMatrix_libMesh(DM dm, Mat * A)
{
  PetscFunctionBegin;
  PetscErrorCode ierr;
  DM_libMesh     * dlm = (DM_libMesh *)(dm->data);
  PetscBool eq;

  ierr = PetscObjectTypeCompare((PetscObject)dm, DMLIBMESH, &eq); CHKERRQ(ierr);

  if (!eq)
    LIBMESH_SETERRQ2(((PetscObject)dm)->comm, PETSC_ERR_ARG_WRONG, "DM of type %s, not of type %s", ((PetscObject)dm)->type_name, DMLIBMESH);

  if (!dlm->sys)
    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONGSTATE, "No libMesh system set for DM_libMesh");

  *A = (dynamic_cast<PetscMatrix<Number> *>(dlm->sys->matrix))->mat();
  ierr = PetscObjectReference((PetscObject)(*A)); CHKERRQ(ierr);
  PetscFunctionReturn(LIBMESH_PETSC_SUCCESS);
}


#undef __FUNCT__
#define __FUNCT__ "DMView_libMesh"
static PetscErrorCode  DMView_libMesh(DM dm, PetscViewer viewer)
{
  PetscErrorCode ierr;
  PetscBool isascii;
  const char * name, * prefix;
  DM_libMesh * dlm = (DM_libMesh *)dm->data;
  PetscFunctionBegin;
  ierr = PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii); CHKERRQ(ierr);
  if (isascii) {
    ierr = PetscObjectGetName((PetscObject)dm, &name);     CHKERRQ(ierr);
    ierr = PetscObjectGetOptionsPrefix((PetscObject)dm, &prefix); CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer, "DM libMesh with name %s and prefix %s\n", name, prefix); CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer, "blocks:"); CHKERRQ(ierr);
    std::map<std::string,unsigned int>::iterator bit = dlm->blockids->begin();
    std::map<std::string,unsigned int>::const_iterator bend = dlm->blockids->end();
    for (; bit != bend; ++bit) {
      ierr = PetscViewerASCIIPrintf(viewer, "(%s,%d) ", bit->first.c_str(), bit->second); CHKERRQ(ierr);
    }
    ierr = PetscViewerASCIIPrintf(viewer, "\n"); CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer, "variables:"); CHKERRQ(ierr);
    std::map<std::string,unsigned int>::iterator vit = dlm->varids->begin();
    std::map<std::string,unsigned int>::const_iterator vend = dlm->varids->end();
    for (; vit != vend; ++vit) {
      ierr = PetscViewerASCIIPrintf(viewer, "(%s,%d) ", vit->first.c_str(), vit->second); CHKERRQ(ierr);
    }
    ierr = PetscViewerASCIIPrintf(viewer, "\n"); CHKERRQ(ierr);
    if (dlm->decomposition_type == DMLIBMESH_NO_DECOMPOSITION) {
      ierr = PetscViewerASCIIPrintf(viewer, "No decomposition\n"); CHKERRQ(ierr);
    }
    else {
      if (dlm->decomposition_type == DMLIBMESH_FIELD_DECOMPOSITION) {
        ierr = PetscViewerASCIIPrintf(viewer, "Field decomposition by variable: "); CHKERRQ(ierr);
      }
      else if (dlm->decomposition_type == DMLIBMESH_DOMAIN_DECOMPOSITION) {
        ierr = PetscViewerASCIIPrintf(viewer, "Domain decomposition by block: "); CHKERRQ(ierr);
      }
      else LIBMESH_SETERRQ1(((PetscObject)dm)->comm, PETSC_ERR_PLIB, "Unexpected decomposition type: %d", dlm->decomposition_type);
      /* FIX: decompositions might have different sizes and components on different ranks. */
      for (auto d : index_range(*dlm->decomposition)) {
        std::set<unsigned int>::iterator dbegin  = (*dlm->decomposition)[d].begin();
        std::set<unsigned int>::iterator dit     = (*dlm->decomposition)[d].begin();
        std::set<unsigned int>::iterator dend    = (*dlm->decomposition)[d].end();
        for (; dit != dend; ++dit) {
          if (dit != dbegin) {
            ierr = PetscViewerASCIIPrintf(viewer, ","); CHKERRQ(ierr);
          }
          ierr = PetscViewerASCIIPrintf(viewer, "%u", *dit); CHKERRQ(ierr);
        }
        ierr = PetscViewerASCIIPrintf(viewer, ";"); CHKERRQ(ierr);
      }
      ierr = PetscViewerASCIIPrintf(viewer, "\n"); CHKERRQ(ierr);
    }
  }

  PetscFunctionReturn(LIBMESH_PETSC_SUCCESS);
}

#undef __FUNCT__
#define __FUNCT__ "DMSetUp_libMesh"
static PetscErrorCode  DMSetUp_libMesh(DM dm)
{
  PetscFunctionBegin;
  PetscErrorCode ierr;
  DM_libMesh     * dlm = (DM_libMesh *)(dm->data);
  PetscBool eq;

  ierr = PetscObjectTypeCompare((PetscObject)dm, DMLIBMESH, &eq); CHKERRQ(ierr);

  if (!eq)
    LIBMESH_SETERRQ2(((PetscObject)dm)->comm, PETSC_ERR_ARG_WRONG, "DM of type %s, not of type %s", ((PetscObject)dm)->type_name, DMLIBMESH);

  if (!dlm->sys)
    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONGSTATE, "No libMesh system set for DM_libMesh");
  /*
    Do not evaluate function, Jacobian or bounds for an embedded DM -- the subproblem might not have enough information for that.
  */
  if (!dlm->embedding) {
    ierr = DMSNESSetFunction(dm, SNESFunction_DMlibMesh, (void *)dm); CHKERRQ(ierr);
    ierr = DMSNESSetJacobian(dm, SNESJacobian_DMlibMesh, (void *)dm); CHKERRQ(ierr);
    if (dlm->sys->nonlinear_solver->bounds || dlm->sys->nonlinear_solver->bounds_object)
      {
        ierr = DMSetVariableBounds(dm, DMVariableBounds_libMesh); CHKERRQ(ierr);
      }
  }
  else {
    /*
      Fow now we don't implement even these, although a linear "Dirichlet" subproblem is well-defined.
      Creating the submatrix, however, might require extracting the submatrix preallocation from an unassembled matrix.
    */
    dm->ops->createglobalvector = 0;
    dm->ops->creatematrix = 0;
  }
  PetscFunctionReturn(LIBMESH_PETSC_SUCCESS);
}




#undef __FUNCT__
#define __FUNCT__ "DMDestroy_libMesh"
static PetscErrorCode  DMDestroy_libMesh(DM dm)
{
  DM_libMesh * dlm = (DM_libMesh *)(dm->data);
  PetscErrorCode ierr;
  PetscFunctionBegin;
  delete dlm->varids;
  delete dlm->varnames;
  delete dlm->blockids;
  delete dlm->blocknames;
  delete dlm->decomposition;
  ierr = ISDestroy(&dlm->embedding); CHKERRQ(ierr);
  ierr = PetscFree(dm->data); CHKERRQ(ierr);

  PetscFunctionReturn(LIBMESH_PETSC_SUCCESS);
}

EXTERN_C_BEGIN
#undef __FUNCT__
#define __FUNCT__ "DMCreate_libMesh"
PetscErrorCode  DMCreate_libMesh(DM dm)
{
  PetscErrorCode ierr;
  DM_libMesh     * dlm;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  ierr = PetscNew(&dlm);CHKERRQ(ierr);
  dm->data = dlm;

  /* DMlibMesh impl */
  dlm->varids     = new(std::map<std::string, unsigned int>);
  dlm->blockids   = new(std::map<std::string, unsigned int>);
  dlm->varnames   = new(std::map<unsigned int, std::string>);
  dlm->blocknames = new(std::map<unsigned int, std::string>);
  dlm->decomposition   = LIBMESH_PETSC_NULLPTR;
  dlm->decomposition_type  = DMLIBMESH_NO_DECOMPOSITION;

  /* DM API */
  dm->ops->createglobalvector = DMCreateGlobalVector_libMesh;
  dm->ops->createlocalvector  = 0; // DMCreateLocalVector_libMesh;
  dm->ops->getcoloring        = 0; // DMGetColoring_libMesh;
  dm->ops->creatematrix       = DMCreateMatrix_libMesh;
  dm->ops->createinterpolation= 0; // DMCreateInterpolation_libMesh;

  dm->ops->refine             = 0; // DMRefine_libMesh;
  dm->ops->coarsen            = 0; // DMCoarsen_libMesh;

  // * dm->ops->getinjection was renamed to dm->ops->createinjection in PETSc 5a84ad338 (5 Jul 2019)
  // * dm->ops-getaggregates was removed in PETSc 97779f9a (5 Jul 2019)
  // * Both changes were merged into PETSc master in 94aad3ce (7 Jul 2019).
#if PETSC_VERSION_LESS_THAN(3,12,0)
  dm->ops->getinjection       = 0; // DMGetInjection_libMesh;
  dm->ops->getaggregates      = 0; // DMGetAggregates_libMesh;
#else
  dm->ops->createinjection = 0;
#endif


  dm->ops->createfielddecomposition    = DMCreateFieldDecomposition_libMesh;
  dm->ops->createdomaindecomposition   = DMCreateDomainDecomposition_libMesh;

  dm->ops->destroy            = DMDestroy_libMesh;
  dm->ops->view               = DMView_libMesh;
  dm->ops->setfromoptions     = 0; // DMSetFromOptions_libMesh;
  dm->ops->setup              = DMSetUp_libMesh;

  /* DMlibMesh API */
  ierr = PetscObjectComposeFunction((PetscObject)dm,"DMlibMeshSetSystem_C",DMlibMeshSetSystem_libMesh);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)dm,"DMlibMeshGetSystem_C",DMlibMeshGetSystem_libMesh);CHKERRQ(ierr);

  PetscFunctionReturn(LIBMESH_PETSC_SUCCESS);
}
EXTERN_C_END

#endif // LIBMESH_HAVE_PETSC
