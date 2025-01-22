// The libMesh Finite Element Library.
// Copyright (C) 2002-2025 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "libmesh/petsc_matrix_base.h"
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

  PetscFunctionBegin;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  PetscBool islibmesh;
  LibmeshPetscCallQ(PetscObjectTypeCompare((PetscObject)dm, DMLIBMESH,&islibmesh));
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
  LibmeshPetscCallQ(DMlibMeshSetUpName_Private(dm));
  PetscFunctionReturn(LIBMESH_PETSC_SUCCESS);
}

#undef  __FUNCT__
#define __FUNCT__ "DMlibMeshGetSystem_libMesh"
PetscErrorCode DMlibMeshGetSystem_libMesh(DM dm, NonlinearImplicitSystem *& sys)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  PetscBool islibmesh;
  LibmeshPetscCallQ(PetscObjectTypeCompare((PetscObject)dm, DMLIBMESH,&islibmesh));
  if (!islibmesh) LIBMESH_SETERRQ2(((PetscObject)dm)->comm, PETSC_ERR_ARG_WRONG, "Got DM of type %s, not of type %s", ((PetscObject)dm)->type_name, DMLIBMESH);
  DM_libMesh * dlm = (DM_libMesh *)(dm->data);
  sys = dlm->sys;
  PetscFunctionReturn(LIBMESH_PETSC_SUCCESS);
}


#undef  __FUNCT__
#define __FUNCT__ "DMlibMeshGetBlocks"
PetscErrorCode DMlibMeshGetBlocks(DM dm, PetscInt * n, char *** blocknames)
{
  PetscInt i;
  PetscFunctionBegin;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  PetscBool islibmesh;
  LibmeshPetscCallQ(PetscObjectTypeCompare((PetscObject)dm, DMLIBMESH,&islibmesh));
  if (!islibmesh) LIBMESH_SETERRQ2(((PetscObject)dm)->comm, PETSC_ERR_ARG_WRONG, "Got DM of type %s, not of type %s", ((PetscObject)dm)->type_name, DMLIBMESH);
  DM_libMesh * dlm = (DM_libMesh *)(dm->data);
#if PETSC_RELEASE_GREATER_EQUALS(3,20,0)
  PetscAssertPointer(n,2);
#else
  PetscValidPointer(n,2);
#endif
  *n = cast_int<unsigned int>(dlm->blockids->size());
  if (!blocknames) PetscFunctionReturn(LIBMESH_PETSC_SUCCESS);
  LibmeshPetscCallQ(PetscMalloc(*n*sizeof(char *), blocknames));
  i = 0;
  for (const auto & pr : *(dlm->blockids))
    {
      LibmeshPetscCallQ(PetscStrallocpy(pr.first.c_str(), *blocknames+i));
      ++i;
    }
  PetscFunctionReturn(LIBMESH_PETSC_SUCCESS);
}

#undef  __FUNCT__
#define __FUNCT__ "DMlibMeshGetVariables"
PetscErrorCode DMlibMeshGetVariables(DM dm, PetscInt * n, char *** varnames)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  PetscBool islibmesh;
  PetscInt i;
  LibmeshPetscCallQ(PetscObjectTypeCompare((PetscObject)dm, DMLIBMESH,&islibmesh));
  if (!islibmesh) LIBMESH_SETERRQ2(((PetscObject)dm)->comm, PETSC_ERR_ARG_WRONG, "Got DM of type %s, not of type %s", ((PetscObject)dm)->type_name, DMLIBMESH);
  DM_libMesh * dlm = (DM_libMesh *)(dm->data);
#if PETSC_RELEASE_GREATER_EQUALS(3,20,0)
  PetscAssertPointer(n,2);
#else
  PetscValidPointer(n,2);
#endif
  *n = cast_int<unsigned int>(dlm->varids->size());
  if (!varnames) PetscFunctionReturn(LIBMESH_PETSC_SUCCESS);
  LibmeshPetscCallQ(PetscMalloc(*n*sizeof(char *), varnames));
  i = 0;
  for (const auto & pr : *(dlm->varids))
    {
      LibmeshPetscCallQ(PetscStrallocpy(pr.first.c_str(), *varnames+i));
      ++i;
    }
  PetscFunctionReturn(LIBMESH_PETSC_SUCCESS);
}

#undef  __FUNCT__
#define __FUNCT__ "DMlibMeshSetUpName_Private"
PetscErrorCode DMlibMeshSetUpName_Private(DM dm)
{
  DM_libMesh * dlm = (DM_libMesh *)dm->data;
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
  LibmeshPetscCallQ(PetscObjectSetName((PetscObject)dm, name.c_str()));
  PetscFunctionReturn(LIBMESH_PETSC_SUCCESS);
}


#undef __FUNCT__
#define __FUNCT__ "DMCreateFieldDecomposition_libMesh"
static PetscErrorCode  DMCreateFieldDecomposition_libMesh(DM dm, PetscInt * len, char *** namelist, IS ** islist, DM ** dmlist)
{
  PetscFunctionBegin;
  DM_libMesh     * dlm = (DM_libMesh *)(dm->data);
  NonlinearImplicitSystem * sys = dlm->sys;
  IS emb;
  if (dlm->decomposition_type != DMLIBMESH_FIELD_DECOMPOSITION) PetscFunctionReturn(LIBMESH_PETSC_SUCCESS);

  *len = cast_int<unsigned int>(dlm->decomposition->size());
  if (namelist) {LibmeshPetscCallQ(PetscMalloc(*len*sizeof(char *), namelist));}
  if (islist)   {LibmeshPetscCallQ(PetscMalloc(*len*sizeof(IS),    islist));}
  if (dmlist)   {LibmeshPetscCallQ(PetscMalloc(*len*sizeof(DM),    dmlist));}
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
      LibmeshPetscCallQ(PetscStrallocpy(dname.c_str(),(*namelist)+d));
    }
    if (islist) {
      IS dis;
      PetscInt * darray;
      LibmeshPetscCallQ(PetscMalloc(sizeof(PetscInt)*dindices.size(), &darray));
      numeric_index_type i = 0;
      for (const auto & id : dindices) {
        darray[i] = id;
        ++i;
      }
      LibmeshPetscCallQ(ISCreateGeneral(((PetscObject)dm)->comm,
                                        cast_int<PetscInt>(dindices.size()),
                                        darray, PETSC_OWN_POINTER, &dis));
      if (dlm->embedding) {
        /* Create a relative embedding into the parent's index space. */
        LibmeshPetscCallQ(ISEmbed(dis,dlm->embedding, PETSC_TRUE, &emb));
        PetscInt elen, dlen;
        LibmeshPetscCallQ(ISGetLocalSize(emb, &elen));
        LibmeshPetscCallQ(ISGetLocalSize(dis, &dlen));
        if (elen != dlen) LIBMESH_SETERRQ1(((PetscObject)dm)->comm, PETSC_ERR_PLIB, "Failed to embed subdomain %zu", d);
        LibmeshPetscCallQ(ISDestroy(&dis));
        dis = emb;
      }
      else {
        emb = dis;
      }
      (*islist)[d] = dis;
    }
    if (dmlist) {
      DM ddm;
      LibmeshPetscCallQ(DMCreate(((PetscObject)dm)->comm, &ddm));
      LibmeshPetscCallQ(DMSetType(ddm, DMLIBMESH));
      DM_libMesh * ddlm = (DM_libMesh *)(ddm->data);
      ddlm->sys = dlm->sys;
      /* copy over the block ids and names */
      *ddlm->blockids = *dlm->blockids;
      *ddlm->blocknames = *dlm->blocknames;
      /* set the vars from the d-th part of the decomposition. */
      *ddlm->varids     = dvarids;
      *ddlm->varnames   = dvarnames;
      LibmeshPetscCallQ(PetscObjectReference((PetscObject)emb));
      ddlm->embedding = emb;
      ddlm->embedding_type = DMLIBMESH_FIELD_EMBEDDING;

      LibmeshPetscCallQ(DMlibMeshSetUpName_Private(ddm));
      LibmeshPetscCallQ(DMSetFromOptions(ddm));
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
  DM_libMesh     * dlm = (DM_libMesh *)(dm->data);
  NonlinearImplicitSystem * sys = dlm->sys;
  IS emb;
  if (dlm->decomposition_type != DMLIBMESH_DOMAIN_DECOMPOSITION) PetscFunctionReturn(LIBMESH_PETSC_SUCCESS);
  *len = cast_int<unsigned int>(dlm->decomposition->size());
  if (namelist)      {LibmeshPetscCallQ(PetscMalloc(*len*sizeof(char *), namelist));}
  if (innerislist)   {LibmeshPetscCallQ(PetscMalloc(*len*sizeof(IS),    innerislist));}
  if (outerislist)   *outerislist = LIBMESH_PETSC_NULLPTR; /* FIX: allow mesh-based overlap. */
  if (dmlist)        {LibmeshPetscCallQ(PetscMalloc(*len*sizeof(DM),    dmlist));}
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
      LibmeshPetscCallQ(PetscStrallocpy(dname.c_str(),(*namelist)+d));
    }
    if (innerislist) {
      PetscInt * darray;
      IS dis;
      LibmeshPetscCallQ(PetscMalloc(sizeof(PetscInt)*dindices.size(), &darray));
      numeric_index_type i = 0;
      for (const auto & id : dindices) {
        darray[i] = id;
        ++i;
      }
      LibmeshPetscCallQ(ISCreateGeneral(((PetscObject)dm)->comm,
                                        cast_int<PetscInt>(dindices.size()),
                                        darray, PETSC_OWN_POINTER, &dis));
      if (dlm->embedding) {
        /* Create a relative embedding into the parent's index space. */
        LibmeshPetscCallQ(ISEmbed(dis,dlm->embedding, PETSC_TRUE, &emb));
        PetscInt elen, dlen;
        LibmeshPetscCallQ(ISGetLocalSize(emb, &elen));
        LibmeshPetscCallQ(ISGetLocalSize(dis, &dlen));
        if (elen != dlen) LIBMESH_SETERRQ1(((PetscObject)dm)->comm, PETSC_ERR_PLIB, "Failed to embed field %zu" , d);
        LibmeshPetscCallQ(ISDestroy(&dis));
        dis = emb;
      }
      else {
        emb = dis;
      }
      if (innerislist) {
        LibmeshPetscCallQ(PetscObjectReference((PetscObject)dis));
        (*innerislist)[d] = dis;
      }
      LibmeshPetscCallQ(ISDestroy(&dis));
    }
    if (dmlist) {
      DM ddm;
      LibmeshPetscCallQ(DMCreate(((PetscObject)dm)->comm, &ddm));
      LibmeshPetscCallQ(DMSetType(ddm, DMLIBMESH));
      DM_libMesh * ddlm = (DM_libMesh *)(ddm->data);
      ddlm->sys = dlm->sys;
      /* copy over the varids and varnames */
      *ddlm->varids    = *dlm->varids;
      *ddlm->varnames  = *dlm->varnames;
      /* set the blocks from the d-th part of the decomposition. */
      *ddlm->blockids    = dblockids;
      *ddlm->blocknames  = dblocknames;
      LibmeshPetscCallQ(PetscObjectReference((PetscObject)emb));
      ddlm->embedding = emb;
      ddlm->embedding_type = DMLIBMESH_DOMAIN_EMBEDDING;

      LibmeshPetscCallQ(DMlibMeshSetUpName_Private(ddm));
      LibmeshPetscCallQ(DMSetFromOptions(ddm));
      (*dmlist)[d] = ddm;
    }
  }
  PetscFunctionReturn(LIBMESH_PETSC_SUCCESS);
}


#undef  __FUNCT__
#define __FUNCT__ "DMlibMeshCreateFieldDecompositionDM"
PetscErrorCode DMlibMeshCreateFieldDecompositionDM(DM dm, PetscInt dnumber, PetscInt * dsizes, char *** dvarlists, DM * ddm)
{
  PetscBool islibmesh;
  PetscFunctionBegin;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  LibmeshPetscCallQ(PetscObjectTypeCompare((PetscObject)dm, DMLIBMESH,&islibmesh));
  if (!islibmesh) LIBMESH_SETERRQ2(((PetscObject)dm)->comm, PETSC_ERR_ARG_WRONG, "Got DM of type %s, not of type %s", ((PetscObject)dm)->type_name, DMLIBMESH);
  if (dnumber < 0) LIBMESH_SETERRQ1(((PetscObject)dm)->comm, PETSC_ERR_ARG_WRONG, "Negative number %" LIBMESH_PETSCINT_FMT " of decomposition parts", dnumber);
#if PETSC_RELEASE_GREATER_EQUALS(3,20,0)
  PetscAssertPointer(ddm,5);
#else
  PetscValidPointer(ddm,5);
#endif
  DM_libMesh * dlm = (DM_libMesh *)(dm->data);
  LibmeshPetscCallQ(DMCreate(((PetscObject)dm)->comm, ddm));
  LibmeshPetscCallQ(DMSetType(*ddm, DMLIBMESH));
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
  LibmeshPetscCallQ(DMlibMeshSetUpName_Private(*ddm));
  LibmeshPetscCallQ(DMSetFromOptions(*ddm));
  LibmeshPetscCallQ(DMSetUp(*ddm));
  PetscFunctionReturn(LIBMESH_PETSC_SUCCESS);
}

#undef  __FUNCT__
#define __FUNCT__ "DMlibMeshCreateDomainDecompositionDM"
PetscErrorCode DMlibMeshCreateDomainDecompositionDM(DM dm, PetscInt dnumber, PetscInt * dsizes, char *** dblocklists, DM * ddm)
{
  PetscBool islibmesh;
  PetscFunctionBegin;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  LibmeshPetscCallQ(PetscObjectTypeCompare((PetscObject)dm, DMLIBMESH,&islibmesh));
  if (!islibmesh) LIBMESH_SETERRQ2(((PetscObject)dm)->comm, PETSC_ERR_ARG_WRONG, "Got DM of type %s, not of type %s", ((PetscObject)dm)->type_name, DMLIBMESH);
  if (dnumber < 0) LIBMESH_SETERRQ1(((PetscObject)dm)->comm, PETSC_ERR_ARG_WRONG, "Negative number %" LIBMESH_PETSCINT_FMT " of decomposition parts", dnumber);
#if PETSC_RELEASE_GREATER_EQUALS(3,20,0)
  PetscAssertPointer(ddm,5);
#else
  PetscValidPointer(ddm,5);
#endif
  DM_libMesh * dlm = (DM_libMesh *)(dm->data);
  LibmeshPetscCallQ(DMCreate(((PetscObject)dm)->comm, ddm));
  LibmeshPetscCallQ(DMSetType(*ddm, DMLIBMESH));
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
  LibmeshPetscCallQ(DMlibMeshSetUpName_Private(*ddm));
  LibmeshPetscCallQ(DMSetFromOptions(*ddm));
  LibmeshPetscCallQ(DMSetUp(*ddm));
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
  PetscFunctionBegin;
  libmesh_assert(x);
  libmesh_assert(r);

  NonlinearImplicitSystem * _sys;
  LibmeshPetscCallQ(DMlibMeshGetSystem(dm, _sys));
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
  PetscFunctionBegin;
  LibmeshPetscCallQ(DMlibMeshFunction(dm,x,r));
  PetscFunctionReturn(LIBMESH_PETSC_SUCCESS);
}


#undef __FUNCT__
#define __FUNCT__ "DMlibMeshJacobian"
static PetscErrorCode DMlibMeshJacobian(DM dm, Vec x, Mat jac, Mat pc)
{
  PetscFunctionBegin;
  NonlinearImplicitSystem * _sys;
  LibmeshPetscCallQ(DMlibMeshGetSystem(dm, _sys));
  NonlinearImplicitSystem & sys = *_sys;

  libmesh_assert(pc);
  libmesh_assert(jac);
  PetscMatrixBase<Number> & the_pc = *PetscMatrixBase<Number>::get_context(pc, sys.comm());
  PetscMatrixBase<Number> & Jac = *PetscMatrixBase<Number>::get_context(jac, sys.comm());
  PetscVector<Number> & X_sys = *cast_ptr<PetscVector<Number> *>(sys.solution.get());
  PetscMatrixBase<Number> & Jac_sys = *cast_ptr<PetscMatrixBase<Number> *>(sys.matrix);
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
  PetscFunctionBegin;
  LibmeshPetscCallQ(DMlibMeshJacobian(dm,x,jac,pc));
  PetscFunctionReturn(LIBMESH_PETSC_SUCCESS);
}

#undef __FUNCT__
#define __FUNCT__ "DMVariableBounds_libMesh"
static PetscErrorCode DMVariableBounds_libMesh(DM dm, Vec xl, Vec xu)
{
  NonlinearImplicitSystem * _sys;
  LibmeshPetscCallQ(DMlibMeshGetSystem(dm, _sys));
  NonlinearImplicitSystem & sys = *_sys;
  PetscVector<Number> XL(xl, sys.comm());
  PetscVector<Number> XU(xu, sys.comm());
  PetscFunctionBegin;
  // Workaround for nonstandard Q suffix warning with quad precision
  const PetscReal petsc_inf = std::numeric_limits<PetscReal>::max() / 4;
  LibmeshPetscCallQ(VecSet(xl, -petsc_inf));
  LibmeshPetscCallQ(VecSet(xu, petsc_inf));
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
  DM_libMesh     * dlm = (DM_libMesh *)(dm->data);
  PetscBool eq;

  LibmeshPetscCallQ(PetscObjectTypeCompare((PetscObject)dm, DMLIBMESH, &eq));

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
    LibmeshPetscCallQ(VecCreate(((PetscObject)v)->comm, x));
    LibmeshPetscCallQ(ISGetLocalSize(dlm->embedding, &n));
    LibmeshPetscCallQ(VecSetSizes(*x,n,PETSC_DETERMINE));
    LibmeshPetscCallQ(VecSetType(*x,((PetscObject)v)->type_name));
    LibmeshPetscCallQ(VecSetFromOptions(*x));
    LibmeshPetscCallQ(VecSetUp(*x));
  }
  else {
    LibmeshPetscCallQ(VecDuplicate(v,x));
  }

#if PETSC_VERSION_LESS_THAN(3,13,0)
  LibmeshPetscCallQ(PetscObjectCompose((PetscObject)*x,"DM",(PetscObject)dm));
#else
  LibmeshPetscCallQ(VecSetDM(*x, dm));
#endif

  PetscFunctionReturn(LIBMESH_PETSC_SUCCESS);
}




#undef __FUNCT__
#define __FUNCT__ "DMCreateMatrix_libMesh"
static PetscErrorCode DMCreateMatrix_libMesh(DM dm, Mat * A)
{
  PetscFunctionBegin;
  DM_libMesh     * dlm = (DM_libMesh *)(dm->data);
  PetscBool eq;

  LibmeshPetscCallQ(PetscObjectTypeCompare((PetscObject)dm, DMLIBMESH, &eq));

  if (!eq)
    LIBMESH_SETERRQ2(((PetscObject)dm)->comm, PETSC_ERR_ARG_WRONG, "DM of type %s, not of type %s", ((PetscObject)dm)->type_name, DMLIBMESH);

  if (!dlm->sys)
    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONGSTATE, "No libMesh system set for DM_libMesh");

  *A = (dynamic_cast<PetscMatrixBase<Number> *>(dlm->sys->matrix))->mat();
  LibmeshPetscCallQ(PetscObjectReference((PetscObject)(*A)));
  PetscFunctionReturn(LIBMESH_PETSC_SUCCESS);
}


#undef __FUNCT__
#define __FUNCT__ "DMView_libMesh"
static PetscErrorCode  DMView_libMesh(DM dm, PetscViewer viewer)
{
  PetscBool isascii;
  const char * name, * prefix;
  DM_libMesh * dlm = (DM_libMesh *)dm->data;
  PetscFunctionBegin;
  LibmeshPetscCallQ(PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii));
  if (isascii) {
    LibmeshPetscCallQ(PetscObjectGetName((PetscObject)dm, &name));
    LibmeshPetscCallQ(PetscObjectGetOptionsPrefix((PetscObject)dm, &prefix));
    LibmeshPetscCallQ(PetscViewerASCIIPrintf(viewer, "DM libMesh with name %s and prefix %s\n", name, prefix));
    LibmeshPetscCallQ(PetscViewerASCIIPrintf(viewer, "blocks:"));
    std::map<std::string,unsigned int>::iterator bit = dlm->blockids->begin();
    std::map<std::string,unsigned int>::const_iterator bend = dlm->blockids->end();
    for (; bit != bend; ++bit) {
      LibmeshPetscCallQ(PetscViewerASCIIPrintf(viewer, "(%s,%d) ", bit->first.c_str(), bit->second));
    }
    LibmeshPetscCallQ(PetscViewerASCIIPrintf(viewer, "\n"));
    LibmeshPetscCallQ(PetscViewerASCIIPrintf(viewer, "variables:"));
    std::map<std::string,unsigned int>::iterator vit = dlm->varids->begin();
    std::map<std::string,unsigned int>::const_iterator vend = dlm->varids->end();
    for (; vit != vend; ++vit) {
      LibmeshPetscCallQ(PetscViewerASCIIPrintf(viewer, "(%s,%d) ", vit->first.c_str(), vit->second));
    }
    LibmeshPetscCallQ(PetscViewerASCIIPrintf(viewer, "\n"));
    if (dlm->decomposition_type == DMLIBMESH_NO_DECOMPOSITION) {
      LibmeshPetscCallQ(PetscViewerASCIIPrintf(viewer, "No decomposition\n"));
    }
    else {
      if (dlm->decomposition_type == DMLIBMESH_FIELD_DECOMPOSITION) {
        LibmeshPetscCallQ(PetscViewerASCIIPrintf(viewer, "Field decomposition by variable: "));
      }
      else if (dlm->decomposition_type == DMLIBMESH_DOMAIN_DECOMPOSITION) {
        LibmeshPetscCallQ(PetscViewerASCIIPrintf(viewer, "Domain decomposition by block: "));
      }
      else LIBMESH_SETERRQ1(((PetscObject)dm)->comm, PETSC_ERR_PLIB, "Unexpected decomposition type: %d", dlm->decomposition_type);
      /* FIX: decompositions might have different sizes and components on different ranks. */
      for (auto d : index_range(*dlm->decomposition)) {
        std::set<unsigned int>::iterator dbegin  = (*dlm->decomposition)[d].begin();
        std::set<unsigned int>::iterator dit     = (*dlm->decomposition)[d].begin();
        std::set<unsigned int>::iterator dend    = (*dlm->decomposition)[d].end();
        for (; dit != dend; ++dit) {
          if (dit != dbegin) {
            LibmeshPetscCallQ(PetscViewerASCIIPrintf(viewer, ","));
          }
          LibmeshPetscCallQ(PetscViewerASCIIPrintf(viewer, "%u", *dit));
        }
        LibmeshPetscCallQ(PetscViewerASCIIPrintf(viewer, ";"));
      }
      LibmeshPetscCallQ(PetscViewerASCIIPrintf(viewer, "\n"));
    }
  }

  PetscFunctionReturn(LIBMESH_PETSC_SUCCESS);
}

#undef __FUNCT__
#define __FUNCT__ "DMSetUp_libMesh"
static PetscErrorCode  DMSetUp_libMesh(DM dm)
{
  PetscFunctionBegin;
  DM_libMesh     * dlm = (DM_libMesh *)(dm->data);
  PetscBool eq;

  LibmeshPetscCallQ(PetscObjectTypeCompare((PetscObject)dm, DMLIBMESH, &eq));

  if (!eq)
    LIBMESH_SETERRQ2(((PetscObject)dm)->comm, PETSC_ERR_ARG_WRONG, "DM of type %s, not of type %s", ((PetscObject)dm)->type_name, DMLIBMESH);

  if (!dlm->sys)
    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONGSTATE, "No libMesh system set for DM_libMesh");
  /*
    Do not evaluate function, Jacobian or bounds for an embedded DM -- the subproblem might not have enough information for that.
  */
  if (!dlm->embedding) {
    LibmeshPetscCallQ(DMSNESSetFunction(dm, SNESFunction_DMlibMesh, (void *)dm));
    LibmeshPetscCallQ(DMSNESSetJacobian(dm, SNESJacobian_DMlibMesh, (void *)dm));
    if (dlm->sys->nonlinear_solver->bounds || dlm->sys->nonlinear_solver->bounds_object)
      LibmeshPetscCallQ(DMSetVariableBounds(dm, DMVariableBounds_libMesh));
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
  PetscFunctionBegin;
  delete dlm->varids;
  delete dlm->varnames;
  delete dlm->blockids;
  delete dlm->blocknames;
  delete dlm->decomposition;
  LibmeshPetscCallQ(ISDestroy(&dlm->embedding));
  LibmeshPetscCallQ(PetscFree(dm->data));

  PetscFunctionReturn(LIBMESH_PETSC_SUCCESS);
}

EXTERN_C_BEGIN
#undef __FUNCT__
#define __FUNCT__ "DMCreate_libMesh"
PetscErrorCode  DMCreate_libMesh(DM dm)
{
  DM_libMesh     * dlm;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  LibmeshPetscCallQ(PetscNew(&dlm));
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
  LibmeshPetscCallQ(PetscObjectComposeFunction((PetscObject)dm,"DMlibMeshSetSystem_C",DMlibMeshSetSystem_libMesh));
  LibmeshPetscCallQ(PetscObjectComposeFunction((PetscObject)dm,"DMlibMeshGetSystem_C",DMlibMeshGetSystem_libMesh));

  PetscFunctionReturn(LIBMESH_PETSC_SUCCESS);
}
EXTERN_C_END

#endif // LIBMESH_HAVE_PETSC
