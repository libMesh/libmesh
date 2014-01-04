#include "libmesh/petsc_macro.h"
// This only works with petsc-3.3 and above.

#if !PETSC_VERSION_LESS_THAN(3,3,0)

// PETSc includes
#include <petsc-private/dmimpl.h>

// Local Includes
#include "libmesh/libmesh_common.h"
#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/petsc_dm_nonlinear_solver.h"
#include "libmesh/petsc_linear_solver.h"
#include "libmesh/petsc_vector.h"
#include "libmesh/petsc_matrix.h"
#include "libmesh/petscdmlibmesh.h"
#include "libmesh/dof_map.h"
#include "libmesh/preconditioner.h"

using namespace libMesh;

#define DMLIBMESH_NO_DECOMPOSITION     0
#define DMLIBMESH_FIELD_DECOMPOSITION  1
#define DMLIBMESH_DOMAIN_DECOMPOSITION 2

#define DMLIBMESH_NO_EMBEDDING         0
#define DMLIBMESH_FIELD_EMBEDDING      1
#define DMLIBMESH_DOMAIN_EMBEDDING     2

struct DM_libMesh
{
  NonlinearImplicitSystem* sys;
  std::map<std::string, unsigned int>  *varids;
  std::map<unsigned int, std::string>  *varnames;
  std::map<std::string, unsigned int>  *blockids;
  std::map<unsigned int, std::string>  *blocknames;
  unsigned int                          decomposition_type;
  std::vector<std::set<unsigned int> > *decomposition;
  unsigned int                          embedding_type;
  IS                                    embedding;
};

#undef  __FUNCT__
#define __FUNCT__ "DMLibMeshGetBlocks"
PetscErrorCode DMLibMeshGetBlocks(DM dm, PetscInt *n, char*** blocknames)
{
  PetscErrorCode ierr;
  PetscInt i;
  PetscFunctionBegin;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  PetscBool islibmesh;
  ierr = PetscObjectTypeCompare((PetscObject)dm, DMLIBMESH,&islibmesh);
  if(!islibmesh) SETERRQ2(((PetscObject)dm)->comm, PETSC_ERR_ARG_WRONG, "Got DM oftype %s, not of type %s", ((PetscObject)dm)->type_name, DMLIBMESH);
  DM_libMesh *dlm = (DM_libMesh *)(dm->data);
  PetscValidPointer(n,2);
  *n = dlm->blockids->size();
  if(!blocknames) PetscFunctionReturn(0);
  ierr = PetscMalloc(*n*sizeof(char*), blocknames); CHKERRQ(ierr);
  i = 0;
  for(std::map<std::string, unsigned int>::const_iterator it = dlm->blockids->begin(); it != dlm->blockids->end(); ++it){
    ierr = PetscStrallocpy(it->first.c_str(), *blocknames+i); CHKERRQ(ierr);
    ++i;
  }
  PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "DMLibMeshGetVariables"
PetscErrorCode DMLibMeshGetVariables(DM dm, PetscInt *n, char*** varnames)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  PetscBool islibmesh;
  PetscInt i;
  ierr = PetscObjectTypeCompare((PetscObject)dm, DMLIBMESH,&islibmesh);
  if(!islibmesh) SETERRQ2(((PetscObject)dm)->comm, PETSC_ERR_ARG_WRONG, "Got DM oftype %s, not of type %s", ((PetscObject)dm)->type_name, DMLIBMESH);
  DM_libMesh *dlm = (DM_libMesh *)(dm->data);
  PetscValidPointer(n,2);
  *n = dlm->varids->size();
  if(!varnames) PetscFunctionReturn(0);
  ierr = PetscMalloc(*n*sizeof(char*), varnames); CHKERRQ(ierr);
  i = 0;
  for(std::map<std::string, unsigned int>::const_iterator it = dlm->varids->begin(); it != dlm->varids->end(); ++it){
    ierr = PetscStrallocpy(it->first.c_str(), *varnames+i); CHKERRQ(ierr);
    ++i;
  }
  PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "DMLibMeshSetUpName_Private"
PetscErrorCode DMLibMeshSetUpName_Private(DM dm)
{
  DM_libMesh* dlm = (DM_libMesh*)dm->data;
  PetscErrorCode ierr;
  PetscFunctionBegin;
  std::string name = dlm->sys->name();
  std::map<unsigned int, std::string> *dnames = PETSC_NULL, *enames = PETSC_NULL;
  if(dlm->decomposition_type == DMLIBMESH_FIELD_DECOMPOSITION) {
    name += ":dec:var:";
    dnames = dlm->varnames;
  }
  if(dlm->decomposition_type == DMLIBMESH_DOMAIN_DECOMPOSITION) {
    name += ":dec:block:";
    dnames = dlm->blocknames;
  }
  if(dnames) {
    for(unsigned int d = 0; d < dlm->decomposition->size(); ++d) {
      for(std::set<unsigned int>::iterator dit = (*dlm->decomposition)[d].begin(); dit != (*dlm->decomposition)[d].end(); ++dit) {
	unsigned int id = *dit;
	if(dit != (*dlm->decomposition)[d].begin())
	  name += ",";
	name += (*dnames)[id];
      }
      name += ";";
    }
  }
  if(dlm->embedding_type == DMLIBMESH_FIELD_EMBEDDING) {
    name += ":emb:var:";
    enames = dlm->varnames;
  }
  if(dlm->embedding_type == DMLIBMESH_DOMAIN_EMBEDDING) {
    name += ":emb:block:";
    enames = dlm->blocknames;
  }
  if(enames) {
    for(std::map<unsigned int, std::string>::iterator eit = enames->begin(); eit != enames->end(); ++eit) {
      std::string ename = eit->second;
      if(eit != enames->begin())
	name += ",";
      name += ename;
    }
    name += ";";
  }
  ierr = PetscObjectSetName((PetscObject)dm, name.c_str()); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "DMLibMeshSetSystem"
PetscErrorCode DMLibMeshSetSystem(DM dm, NonlinearImplicitSystem& sys)
{
  const Parallel::Communicator &comm(sys.comm());

  PetscErrorCode ierr;
  PetscFunctionBegin;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  PetscBool islibmesh;
  ierr = PetscObjectTypeCompare((PetscObject)dm, DMLIBMESH,&islibmesh);
  if(!islibmesh) SETERRQ2(((PetscObject)dm)->comm, PETSC_ERR_ARG_WRONG, "Got DM oftype %s, not of type %s", ((PetscObject)dm)->type_name, DMLIBMESH);

  if(dm->setupcalled) SETERRQ(((PetscObject)dm)->comm, PETSC_ERR_ARG_WRONGSTATE, "Cannot reset the libMesh system after DM has been set up.");
  DM_libMesh *dlm = (DM_libMesh *)(dm->data);
  dlm->sys =&sys;
  /* Initially populate the sets of active blockids and varids using all of the
     existing blocks/variables (only variables are supported at the moment). */
  DofMap& dofmap = dlm->sys->get_dof_map();
  dlm->varids->clear();
  dlm->varnames->clear();
  for(unsigned int v = 0; v < dofmap.n_variables(); ++v) {
    std::string vname = dofmap.variable(v).name();
    dlm->varids->insert(std::pair<std::string,unsigned int>(vname,v));
    dlm->varnames->insert(std::pair<unsigned int,std::string>(v,vname));
  }
  const MeshBase& mesh = dlm->sys->get_mesh();
  dlm->blockids->clear();
  dlm->blocknames->clear();
  std::set<subdomain_id_type> blocks;
  /* The following effectively is a verbatim copy of MeshBase::n_subdomains(). */
  // This requires an inspection on every processor
  libmesh_parallel_only(mesh.comm());
  MeshBase::const_element_iterator       el  = mesh.active_elements_begin();
  const MeshBase::const_element_iterator end = mesh.active_elements_end();
  for (; el!=end; ++el)
    blocks.insert((*el)->subdomain_id());
  // Some subdomains may only live on other processors
  comm.set_union(blocks);

  std::set<subdomain_id_type>::iterator bit = blocks.begin();
  std::set<subdomain_id_type>::iterator bend = blocks.end();
  if(bit == bend) SETERRQ(((PetscObject)dm)->comm, PETSC_ERR_PLIB, "No mesh blocks found.");

  for(; bit != bend; ++bit) {
    subdomain_id_type bid = *bit;
    std::string bname = mesh.subdomain_name(bid);
    if(!bname.length()) {
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
  ierr = DMLibMeshSetUpName_Private(dm); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "DMLibMeshGetSystem"
PetscErrorCode DMLibMeshGetSystem(DM dm, NonlinearImplicitSystem*& sys)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  PetscBool islibmesh;
  ierr = PetscObjectTypeCompare((PetscObject)dm, DMLIBMESH,&islibmesh); CHKERRQ(ierr);
  if(!islibmesh) SETERRQ2(((PetscObject)dm)->comm, PETSC_ERR_ARG_WRONG, "Got DM oftype %s, not of type %s", ((PetscObject)dm)->type_name, DMLIBMESH);
  DM_libMesh *dlm = (DM_libMesh *)(dm->data);
  sys = dlm->sys;
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "DMCreateFieldDecomposition_libMesh"
static PetscErrorCode  DMCreateFieldDecomposition_libMesh(DM dm, PetscInt *len, char ***namelist, IS **islist, DM **dmlist)
{
  PetscFunctionBegin;
  PetscErrorCode ierr;
  DM_libMesh     *dlm = (DM_libMesh *)(dm->data);
  NonlinearImplicitSystem *sys = dlm->sys;
  IS emb;
  if(dlm->decomposition_type != DMLIBMESH_FIELD_DECOMPOSITION) PetscFunctionReturn(0);

  *len = dlm->decomposition->size();
  if(namelist) {ierr = PetscMalloc(*len*sizeof(char*), namelist);  CHKERRQ(ierr);}
  if(islist)   {ierr = PetscMalloc(*len*sizeof(IS),    islist);    CHKERRQ(ierr);}
  if(dmlist)   {ierr = PetscMalloc(*len*sizeof(DM),    dmlist);    CHKERRQ(ierr);}
  DofMap& dofmap = dlm->sys->get_dof_map();
  for(unsigned int d = 0; d < dlm->decomposition->size(); ++d) {
    std::set<numeric_index_type>         dindices;
    std::string                          dname;
    std::map<std::string, unsigned int>  dvarids;
    std::map<unsigned int, std::string>  dvarnames;
    unsigned int                         dvcount = 0;
    for(std::set<unsigned int>::const_iterator dvit = (*dlm->decomposition)[d].begin(); dvit != (*dlm->decomposition)[d].end(); ++dvit){
      unsigned int v = *dvit;
      std::string vname = (*dlm->varnames)[v];
      dvarids.insert(std::pair<std::string, unsigned int>(vname,v));
      dvarnames.insert(std::pair<unsigned int,std::string>(v,vname));
      if(!dvcount) dname = vname;
      else   dname += "_" + vname;
      ++dvcount;
      if(!islist) continue;
      /* Iterate only over this DM's blocks. */
      for(std::map<std::string, unsigned int>::const_iterator bit = dlm->blockids->begin(); bit != dlm->blockids->end(); ++bit) {
	unsigned int b = bit->second;
	MeshBase::const_element_iterator el     = sys->get_mesh().active_local_subdomain_elements_begin(b);
	MeshBase::const_element_iterator end_el = sys->get_mesh().active_local_subdomain_elements_end(b);
	for ( ; el != end_el; ++el) {
	  const Elem* elem = *el;
	  //unsigned int e_subdomain = elem->subdomain_id();
	  std::vector<numeric_index_type> evindices;
	  // Get the degree of freedom indices for the given variable off the current element.
	  dofmap.dof_indices(elem, evindices, v);
	  for(unsigned int i = 0; i < evindices.size(); ++i) {
	    numeric_index_type dof = evindices[i];
	    if(dof >= dofmap.first_dof() && dof < dofmap.end_dof()) /* might want to use variable_first/last_local_dof instead */
	      dindices.insert(dof);
	  }
	}
      }
    }
    if(namelist) {
      ierr = PetscStrallocpy(dname.c_str(),(*namelist)+d);            CHKERRQ(ierr);
    }
    if(islist) {
      IS dis;
      PetscInt *darray;
      ierr = PetscMalloc(sizeof(PetscInt)*dindices.size(), &darray); CHKERRQ(ierr);
      numeric_index_type i = 0;
      for(std::set<numeric_index_type>::const_iterator it = dindices.begin(); it != dindices.end(); ++it) {
	darray[i] = *it;
	++i;
      }
      ierr = ISCreateGeneral(((PetscObject)dm)->comm, dindices.size(),darray, PETSC_OWN_POINTER, &dis); CHKERRQ(ierr);
      if(dlm->embedding) {
	/* Create a relative embedding into the parent's index space. */
#if PETSC_RELEASE_LESS_THAN(3,3,1)
	ierr = ISMapFactorRight(dis,dlm->embedding, PETSC_TRUE, &emb); CHKERRQ(ierr);
#else
	ierr = ISEmbed(dis,dlm->embedding, PETSC_TRUE, &emb); CHKERRQ(ierr);
#endif
	PetscInt elen, dlen;
	ierr = ISGetLocalSize(emb, &elen); CHKERRQ(ierr);
	ierr = ISGetLocalSize(dis, &dlen); CHKERRQ(ierr);
	if(elen != dlen) SETERRQ1(((PetscObject)dm)->comm, PETSC_ERR_PLIB, "Failed to embed subdomain %D", d);
	ierr = ISDestroy(&dis); CHKERRQ(ierr);
	dis = emb;
      }
      else {
	emb = dis;
      }
      (*islist)[d] = dis;
    }
    if(dmlist) {
      DM ddm;
      ierr = DMCreate(((PetscObject)dm)->comm, &ddm); CHKERRQ(ierr);
      ierr = DMSetType(ddm, DMLIBMESH);               CHKERRQ(ierr);
      DM_libMesh *ddlm = (DM_libMesh*)(ddm->data);
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

      ierr = DMLibMeshSetUpName_Private(ddm); CHKERRQ(ierr);
      ierr = DMSetFromOptions(ddm);           CHKERRQ(ierr);
      (*dmlist)[d] = ddm;
    }
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMCreateDomainDecomposition_libMesh"
static PetscErrorCode  DMCreateDomainDecomposition_libMesh(DM dm, PetscInt *len, char ***namelist, IS **innerislist, IS **outerislist, DM **dmlist)
{
  PetscFunctionBegin;
  PetscErrorCode ierr;
  DM_libMesh     *dlm = (DM_libMesh *)(dm->data);
  NonlinearImplicitSystem *sys = dlm->sys;
  IS emb;
  if(dlm->decomposition_type != DMLIBMESH_DOMAIN_DECOMPOSITION) PetscFunctionReturn(0);
  *len = dlm->decomposition->size();
  if(namelist)      {ierr = PetscMalloc(*len*sizeof(char*), namelist);  CHKERRQ(ierr);}
  if(innerislist)   {ierr = PetscMalloc(*len*sizeof(IS),    innerislist);    CHKERRQ(ierr);}
  if(outerislist)   *outerislist = PETSC_NULL; /* FIX: allow mesh-based overlap. */
  if(dmlist)        {ierr = PetscMalloc(*len*sizeof(DM),    dmlist);    CHKERRQ(ierr);}
  for(unsigned int d = 0; d < dlm->decomposition->size(); ++d) {
    std::set<numeric_index_type>               dindices;
    std::string                          dname;
    std::map<std::string, unsigned int>  dblockids;
    std::map<unsigned int,std::string>   dblocknames;
    unsigned int                         dbcount = 0;
    for(std::set<unsigned int>::const_iterator bit = (*dlm->decomposition)[d].begin(); bit != (*dlm->decomposition)[d].end(); ++bit){
      unsigned int b = *bit;
      std::string bname = (*dlm->blocknames)[b];
      dblockids.insert(std::pair<std::string, unsigned int>(bname,b));
      dblocknames.insert(std::pair<unsigned int,std::string>(b,bname));
      if(!dbcount) dname = bname;
      else   dname += "_" + bname;
      ++dbcount;
      if(!innerislist) continue;
      MeshBase::const_element_iterator       el     = sys->get_mesh().active_local_subdomain_elements_begin(b);
      const MeshBase::const_element_iterator end_el = sys->get_mesh().active_local_subdomain_elements_end(b);
      for ( ; el != end_el; ++el) {
	const Elem* elem = *el;
	std::vector<numeric_index_type> evindices;
	/* Iterate only over this DM's variables. */
	for(std::map<std::string, unsigned int>::const_iterator vit = dlm->varids->begin(); vit != dlm->varids->end(); ++vit) {
	  unsigned int v = vit->second;
	  // Get the degree of freedom indices for the given variable off the current element.
	  sys->get_dof_map().dof_indices(elem, evindices, v);
	  for(unsigned int i = 0; i < evindices.size(); ++i) {
	    numeric_index_type dof = evindices[i];
	    if(dof >= sys->get_dof_map().first_dof() && dof < sys->get_dof_map().end_dof()) /* might want to use variable_first/last_local_dof instead */
	      dindices.insert(dof);
	  }
	}
      }
    }
    if(namelist) {
      ierr = PetscStrallocpy(dname.c_str(),(*namelist)+d);            CHKERRQ(ierr);
    }
    if(innerislist) {
      PetscInt *darray;
      IS dis;
      ierr = PetscMalloc(sizeof(PetscInt)*dindices.size(), &darray); CHKERRQ(ierr);
      numeric_index_type i = 0;
      for(std::set<numeric_index_type>::const_iterator it = dindices.begin(); it != dindices.end(); ++it) {
	darray[i] = *it;
	++i;
      }
      ierr = ISCreateGeneral(((PetscObject)dm)->comm, dindices.size(),darray, PETSC_OWN_POINTER, &dis); CHKERRQ(ierr);
      if(dlm->embedding) {
	/* Create a relative embedding into the parent's index space. */
#if PETSC_RELEASE_LESS_THAN(3,3,1)
	ierr = ISMapFactorRight(dis,dlm->embedding, PETSC_TRUE, &emb); CHKERRQ(ierr);
#else
	ierr = ISEmbed(dis,dlm->embedding, PETSC_TRUE, &emb); CHKERRQ(ierr);
#endif
	PetscInt elen, dlen;
	ierr = ISGetLocalSize(emb, &elen); CHKERRQ(ierr);
	ierr = ISGetLocalSize(dis, &dlen);  CHKERRQ(ierr);
	if(elen != dlen) SETERRQ1(((PetscObject)dm)->comm, PETSC_ERR_PLIB, "Failed to embed field %D", d);
	ierr = ISDestroy(&dis); CHKERRQ(ierr);
	dis = emb;
      }
      else {
	emb = dis;
      }
      if(innerislist) {
	ierr = PetscObjectReference((PetscObject)dis); CHKERRQ(ierr);
	(*innerislist)[d] = dis;
      }
      ierr = ISDestroy(&dis); CHKERRQ(ierr);
    }
    if(dmlist) {
      DM ddm;
      ierr = DMCreate(((PetscObject)dm)->comm, &ddm); CHKERRQ(ierr);
      ierr = DMSetType(ddm, DMLIBMESH);               CHKERRQ(ierr);
      DM_libMesh *ddlm = (DM_libMesh*)(ddm->data);
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

      ierr = DMLibMeshSetUpName_Private(ddm); CHKERRQ(ierr);
      ierr = DMSetFromOptions(ddm);           CHKERRQ(ierr);
      (*dmlist)[d] = ddm;
    }
  }
  PetscFunctionReturn(0);
}


#undef  __FUNCT__
#define __FUNCT__ "DMLibMeshCreateFieldDecompositionDM"
PetscErrorCode DMLibMeshCreateFieldDecompositionDM(DM dm, PetscInt dnumber, PetscInt* dsizes, char*** dvarlists, DM* ddm)
{
  PetscErrorCode ierr;
  PetscBool islibmesh;
  PetscFunctionBegin;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  ierr = PetscObjectTypeCompare((PetscObject)dm, DMLIBMESH,&islibmesh);
  if(!islibmesh) SETERRQ2(((PetscObject)dm)->comm, PETSC_ERR_ARG_WRONG, "Got DM oftype %s, not of type %s", ((PetscObject)dm)->type_name, DMLIBMESH);
  if(dnumber < 0) SETERRQ1(((PetscObject)dm)->comm, PETSC_ERR_ARG_WRONG, "Negative number %D of decomposition parts", dnumber);
  PetscValidPointer(ddm,5);
  DM_libMesh *dlm = (DM_libMesh *)(dm->data);
  ierr = DMCreate(((PetscObject)dm)->comm, ddm); CHKERRQ(ierr);
  ierr = DMSetType(*ddm, DMLIBMESH);             CHKERRQ(ierr);
  DM_libMesh *ddlm = (DM_libMesh *)((*ddm)->data);
  ddlm->sys = dlm->sys;
  ddlm->varids = dlm->varids;
  ddlm->varnames = dlm->varnames;
  ddlm->blockids = dlm->blockids;
  ddlm->blocknames = dlm->blocknames;
  ddlm->decomposition = new(std::vector<std::set<unsigned int> >);
  ddlm->decomposition_type = DMLIBMESH_FIELD_DECOMPOSITION;
  if(dnumber) {
    for(PetscInt d = 0; d < dnumber; ++d) {
      if(dsizes[d] < 0) SETERRQ2(((PetscObject)dm)->comm, PETSC_ERR_ARG_WRONG, "Negative size %D of decomposition part %D", dsizes[d],d);
      ddlm->decomposition->push_back(std::set<unsigned int>());
      for(PetscInt v = 0; v < dsizes[d]; ++v) {
	std::string vname(dvarlists[d][v]);
	std::map<std::string, unsigned int>::const_iterator vit = dlm->varids->find(vname);
	if(vit == dlm->varids->end())
	  SETERRQ3(((PetscObject)dm)->comm, PETSC_ERR_ARG_WRONG, "Variable %D on the %D-th list with name %s is not owned by this DM", v, d, dvarlists[d][v]);
	unsigned int vid = vit->second;
	(*ddlm->decomposition)[d].insert(vid);
      }
    }
  }
  else { /* Empty splits indicate default: split all variables with one per split. */
    PetscInt d = 0;
    for(std::map<std::string, unsigned int>::const_iterator vit = ddlm->varids->begin(); vit != ddlm->varids->end(); ++vit) {
      ddlm->decomposition->push_back(std::set<unsigned int>());
      unsigned int vid = vit->second;
      std::string vname = vit->first;
      (*ddlm->decomposition)[d].insert(vid);
      ++d;
    }
  }
  ierr = DMLibMeshSetUpName_Private(*ddm); CHKERRQ(ierr);
  ierr = DMSetFromOptions(*ddm);           CHKERRQ(ierr);
  ierr = DMSetUp(*ddm);                    CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "DMLibMeshCreateDomainDecompositionDM"
PetscErrorCode DMLibMeshCreateDomainDecompositionDM(DM dm, PetscInt dnumber, PetscInt* dsizes, char*** dblocklists, DM* ddm)
{
  PetscErrorCode ierr;
  PetscBool islibmesh;
  PetscFunctionBegin;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  ierr = PetscObjectTypeCompare((PetscObject)dm, DMLIBMESH,&islibmesh);
  if(!islibmesh) SETERRQ2(((PetscObject)dm)->comm, PETSC_ERR_ARG_WRONG, "Got DM oftype %s, not of type %s", ((PetscObject)dm)->type_name, DMLIBMESH);
  if(dnumber < 0) SETERRQ1(((PetscObject)dm)->comm, PETSC_ERR_ARG_WRONG, "Negative number %D of decomposition parts", dnumber);
  PetscValidPointer(ddm,5);
  DM_libMesh *dlm = (DM_libMesh *)(dm->data);
  ierr = DMCreate(((PetscObject)dm)->comm, ddm); CHKERRQ(ierr);
  ierr = DMSetType(*ddm, DMLIBMESH);             CHKERRQ(ierr);
  DM_libMesh *ddlm = (DM_libMesh *)((*ddm)->data);
  ddlm->sys = dlm->sys;
  ddlm->varids   = dlm->varids;
  ddlm->varnames = dlm->varnames;
  ddlm->blockids   = dlm->blockids;
  ddlm->blocknames = dlm->blocknames;
  ddlm->decomposition = new(std::vector<std::set<unsigned int> >);
  ddlm->decomposition_type = DMLIBMESH_DOMAIN_DECOMPOSITION;
  if(dnumber) {
    for(PetscInt d = 0; d < dnumber; ++d) {
      if(dsizes[d] < 0) SETERRQ2(((PetscObject)dm)->comm, PETSC_ERR_ARG_WRONG, "Negative size %D of decomposition part %D", dsizes[d],d);
      ddlm->decomposition->push_back(std::set<unsigned int>());
      for(PetscInt b = 0; b < dsizes[d]; ++b) {
	std::string bname(dblocklists[d][b]);
	std::map<std::string, unsigned int>::const_iterator bit = dlm->blockids->find(bname);
	if(bit == dlm->blockids->end())
	  SETERRQ3(((PetscObject)dm)->comm, PETSC_ERR_ARG_WRONG, "Block %D on the %D-th list with name %s is not owned by this DM", b, d, dblocklists[d][b]);
	unsigned int bid = bit->second;
	(*ddlm->decomposition)[d].insert(bid);
      }
    }
  }
  else { /* Empty splits indicate default: split all blocks with one per split. */
    PetscInt d = 0;
    for(std::map<std::string, unsigned int>::const_iterator bit = ddlm->blockids->begin(); bit != ddlm->blockids->end(); ++bit) {
      ddlm->decomposition->push_back(std::set<unsigned int>());
      unsigned int bid = bit->second;
      std::string bname = bit->first;
      (*ddlm->decomposition)[d].insert(bid);
      ++d;
    }
  }
  ierr = DMLibMeshSetUpName_Private(*ddm); CHKERRQ(ierr);
  ierr = DMSetFromOptions(*ddm);           CHKERRQ(ierr);
  ierr = DMSetUp(*ddm);                    CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

struct token {
  const  char* s;
  struct token *next;
};
#undef __FUNCT__
#define __FUNCT__ "DMLibMeshParseDecompositionDescriptor_Private"
static PetscErrorCode  DMLibMeshParseDecompositionDescriptor_Private(DM dm, const char* ddesc, PetscInt* dtype, PetscInt* dcount, PetscInt** dsizes, char ****dlists)
{
  PetscFunctionBegin;
  PetscErrorCode ierr;
  PetscBool eq;
  char *s0, *s, *ss;
  struct token *llfirst = PETSC_NULL, *lllast = PETSC_NULL, *tok;
  PetscInt stcount = 0, brcount = 0, d, i;
  size_t len0, count;

  /*
     Parse the decomposition descriptor.
     Decomposition names could be of one of two forms:
     var:v1,v2;v3,v4;v4,v5;
     block:b1,b2;b3,b4;b4,b5;
     resulting in an overlapping decomposition that groups
     variables (v1,v2), (v3,v4), (v4,v5) or
     blocks    (b1,b2), (b3,b4), (b4,b5).
  */
  /* Copy the descriptor so that we can manipulate it in place. */
  ierr = PetscStrallocpy(ddesc,&s0);   CHKERRQ(ierr);
  ierr = PetscStrlen(s0, &len0)  ;     CHKERRQ(ierr);
  ierr = PetscStrstr(s0,":",&ss);      CHKERRQ(ierr);
  if(!ss) {
    ss = s0+len0;
  }
  else {
    *ss = 0;
  }
  ierr = PetscStrcmp(s0,"var",&eq);    CHKERRQ(ierr);
  if(eq) {
    *dtype=DMLIBMESH_FIELD_DECOMPOSITION;
  }
  else {
    ierr = PetscStrcmp(s0,"block",&eq);CHKERRQ(ierr);
    if(!eq)
      SETERRQ1(((PetscObject)dm)->comm, PETSC_ERR_ARG_WRONG, "Could not determine decomposition type from descriptor: %s\n", ddesc); CHKERRQ(ierr);
    *dtype=DMLIBMESH_DOMAIN_DECOMPOSITION;
  }
  ierr = PetscStrlen(s0,&count);       CHKERRQ(ierr);
  while(count < len0) {
    struct token *st, *br;
    ++ss; ++count;
    s=ss;
    while(*ss && *ss != ',' && *ss != ';') {
      ++ss; ++count;
    }
    st = PETSC_NULL; br = PETSC_NULL;
    if(*ss) {
      /*
	 Found a separator, or a break.
	 Add an appropriate token to the list.
	 A token separator ',' produces no token.
      */
      if(*ss == ';') {
	/* Create a break token: a token with a null string. */
#if PETSC_RELEASE_LESS_THAN(3,5,0)
	ierr = PetscNew(struct token,&br);CHKERRQ(ierr);
#else
	ierr = PetscNew(&br);CHKERRQ(ierr);
#endif
      }
      *ss = 0;
      if(s != ss) {
	/* A nonempty string. */
#if PETSC_RELEASE_LESS_THAN(3,5,0)
	ierr = PetscNew(struct token, &st);CHKERRQ(ierr);
#else
	ierr = PetscNew(&st);CHKERRQ(ierr);
#endif
	st->s = s; /* The string will be properly copied below. */
      }
      /* Add the new tokens to the list. */
      if(st) {
	if(!lllast) {
	  llfirst = lllast = st;
	}
	else {
	  lllast->next = st; lllast = st;
	}
      }
      if(br) {
	if(!lllast) {
	  llfirst = lllast = br;
	}
	else {
	  lllast->next = br; lllast = br;
	}
      }
    }
  }
  /* The result of parsing is in the linked list ll. */
  /* Count up the strings and the breaks. */
  tok = llfirst;
  while(tok) {
    if(tok->s)
      ++stcount;
    else
      ++brcount;
    tok = tok->next;
  }
  /* Allocate the space for the output. */
  *dcount = brcount;
  ierr = PetscMalloc(*dcount*sizeof(PetscInt), dsizes); CHKERRQ(ierr);
  ierr = PetscMalloc(*dcount*sizeof(char**),   dlists); CHKERRQ(ierr);
  for(d = 0; d < *dcount; ++d) (*dsizes)[d] = 0;
  tok = llfirst; d = 0;
  while(tok) {
    if(tok->s)
      ++(*dsizes)[d];
    else
      ++d;
    tok = tok->next;
  }
  for(d = 0; d < *dcount; ++d) {
    ierr = PetscMalloc(sizeof(char**)*(*dsizes)[d], (*dlists)+d); CHKERRQ(ierr);
  }
  /* Now copy strings and destroy tokens. */
  tok = llfirst; d = 0; i = 0;
  while(tok) {
    if(tok->s) {
      ierr = PetscStrallocpy(tok->s, (*dlists)[d]+i); CHKERRQ(ierr);
      ++i;
    }
    else {
      ++d;
      i = 0;
    }
    llfirst = tok;
    tok = tok->next;
    ierr = PetscFree(llfirst); CHKERRQ(ierr);
  }
  /* Deallocate workspace. */
  ierr = PetscFree(s0); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMCreateFieldDecompositionDM_libMesh"
static PetscErrorCode  DMCreateFieldDecompositionDM_libMesh(DM dm, const char* ddesc, DM *ddm)
{
  PetscFunctionBegin;
  PetscErrorCode ierr;
  PetscInt dtype, dcount, *dsizes;
  char ***dlists;
  PetscFunctionBegin;
  *ddm = PETSC_NULL;
  ierr = DMLibMeshParseDecompositionDescriptor_Private(dm,ddesc,&dtype,&dcount,&dsizes,&dlists); CHKERRQ(ierr);
  if(dtype == DMLIBMESH_FIELD_DECOMPOSITION){
    ierr = DMLibMeshCreateFieldDecompositionDM(dm,dcount,dsizes,dlists,ddm); CHKERRQ(ierr);
  }
  else SETERRQ1(((PetscObject)dm)->comm, PETSC_ERR_PLIB, "Uexpected unknown decomposition type for field decomposition descriptor %s", ddesc);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMCreateDomainDecompositionDM_libMesh"
static PetscErrorCode  DMCreateDomainDecompositionDM_libMesh(DM dm, const char* ddesc, DM *ddm)
{
  PetscFunctionBegin;
  PetscErrorCode ierr;
  PetscInt dtype, dcount, *dsizes;
  char ***dlists;
  PetscFunctionBegin;
  *ddm = PETSC_NULL;
  ierr = DMLibMeshParseDecompositionDescriptor_Private(dm,ddesc,&dtype,&dcount,&dsizes,&dlists); CHKERRQ(ierr);
  if(dtype == DMLIBMESH_DOMAIN_DECOMPOSITION) {
    ierr = DMLibMeshCreateDomainDecompositionDM(dm,dcount,dsizes,dlists,ddm); CHKERRQ(ierr);
  }
  else SETERRQ1(((PetscObject)dm)->comm, PETSC_ERR_PLIB, "Uexpected unknown decomposition type for domain decomposition descriptor %s", ddesc);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMlibMeshFunction"
static PetscErrorCode DMlibMeshFunction(DM dm, Vec x, Vec r)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  libmesh_assert(x);
  libmesh_assert(r);

  NonlinearImplicitSystem* _sys;
  ierr = DMLibMeshGetSystem(dm, _sys); CHKERRQ(ierr);
  NonlinearImplicitSystem& sys = *_sys;
  PetscVector<Number>& X_sys = *libmesh_cast_ptr<PetscVector<Number>* >(sys.solution.get());
  PetscVector<Number>& R_sys = *libmesh_cast_ptr<PetscVector<Number>* >(sys.rhs);
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
  if (_sys->nonlinear_solver->residual && _sys->nonlinear_solver->residual_object)
    {
      libMesh::err << "ERROR: cannot specifiy both a function and object to compute the Residual!" << std::endl;
      libmesh_error();
    }

  if (_sys->nonlinear_solver->matvec && _sys->nonlinear_solver->residual_and_jacobian_object)
    {
      libMesh::err << "ERROR: cannot specifiy both a function and object to compute the combined Residual & Jacobian!" << std::endl;
      libmesh_error();
    }

  if (_sys->nonlinear_solver->residual != NULL)
    _sys->nonlinear_solver->residual(*(_sys->current_local_solution.get()), R, *_sys);

  else if (_sys->nonlinear_solver->residual_object != NULL)
    _sys->nonlinear_solver->residual_object->residual(*(_sys->current_local_solution.get()), R, *_sys);

  else if (_sys->nonlinear_solver->matvec   != NULL)
    _sys->nonlinear_solver->matvec(*(_sys->current_local_solution.get()), &R, NULL, *_sys);

  else if (_sys->nonlinear_solver->residual_and_jacobian_object != NULL)
    _sys->nonlinear_solver->residual_and_jacobian_object->residual_and_jacobian(*(_sys->current_local_solution.get()), &R, NULL, *_sys);

  else
    libmesh_error();

  R.close();
  X_global.close();
  PetscFunctionReturn(0);
}

#if !PETSC_RELEASE_LESS_THAN(3,3,1)
#undef __FUNCT__
#define __FUNCT__ "SNESFunction_DMlibMesh"
static PetscErrorCode SNESFunction_DMlibMesh(SNES, Vec x, Vec r, void *ctx)
{
  DM dm = (DM)ctx;
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = DMlibMeshFunction(dm,x,r);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
#endif


#undef __FUNCT__
#define __FUNCT__ "DMlibMeshJacobian"
static PetscErrorCode DMlibMeshJacobian(DM dm, Vec x, Mat jac, Mat pc, MatStructure *msflag)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  NonlinearImplicitSystem* _sys;
  ierr = DMLibMeshGetSystem(dm, _sys); CHKERRQ(ierr);
  NonlinearImplicitSystem& sys = *_sys;

  PetscMatrix<Number>  the_pc(pc,sys.comm());
  PetscMatrix<Number>  Jac(jac,sys.comm());
  PetscVector<Number>& X_sys = *libmesh_cast_ptr<PetscVector<Number>*>(sys.solution.get());
  PetscMatrix<Number>& Jac_sys = *libmesh_cast_ptr<PetscMatrix<Number>*>(sys.matrix);
  PetscVector<Number>  X_global(x, sys.comm());

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
  if (sys.nonlinear_solver->jacobian && sys.nonlinear_solver->jacobian_object)
    {
      libMesh::err << "ERROR: cannot specifiy both a function and object to compute the Jacobian!" << std::endl;
      libmesh_error();
    }

  if (sys.nonlinear_solver->matvec && sys.nonlinear_solver->residual_and_jacobian_object)
    {
      libMesh::err << "ERROR: cannot specifiy both a function and object to compute the combined Residual & Jacobian!" << std::endl;
      libmesh_error();
    }

  if (sys.nonlinear_solver->jacobian != NULL)
    sys.nonlinear_solver->jacobian(*(sys.current_local_solution.get()), the_pc, sys);

  else if (sys.nonlinear_solver->jacobian_object != NULL)
    sys.nonlinear_solver->jacobian_object->jacobian(*(sys.current_local_solution.get()), the_pc, sys);

  else if (sys.nonlinear_solver->matvec != NULL)
    sys.nonlinear_solver->matvec(*(sys.current_local_solution.get()), NULL, &the_pc, sys);

  else if (sys.nonlinear_solver->residual_and_jacobian_object != NULL)
    sys.nonlinear_solver->residual_and_jacobian_object->residual_and_jacobian(*(sys.current_local_solution.get()), NULL, &the_pc, sys);

  else
    libmesh_error();

  the_pc.close();
  Jac.close();
  X_global.close();

  *msflag = SAME_NONZERO_PATTERN;
  PetscFunctionReturn(0);
}

#if !PETSC_RELEASE_LESS_THAN(3,3,1)
#undef  __FUNCT__
#define __FUNCT__ "SNESJacobian_DMlibMesh"
static PetscErrorCode SNESJacobian_DMlibMesh(SNES,Vec x,Mat *jac,Mat *pc, MatStructure* flag, void* ctx)
{
  DM dm = (DM)ctx;
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = DMlibMeshJacobian(dm,x,*jac,*pc,flag); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
#endif

#undef __FUNCT__
#define __FUNCT__ "DMVariableBounds_libMesh"
static PetscErrorCode DMVariableBounds_libMesh(DM dm, Vec xl, Vec xu)
{
  PetscErrorCode ierr;
  NonlinearImplicitSystem* _sys;
  ierr = DMLibMeshGetSystem(dm, _sys); CHKERRQ(ierr);
  NonlinearImplicitSystem& sys = *_sys;
  PetscVector<Number> XL(xl, sys.comm());
  PetscVector<Number> XU(xu, sys.comm());
  PetscFunctionBegin;

  ierr = VecSet(xl, SNES_VI_NINF); CHKERRQ(ierr);
  ierr = VecSet(xu, SNES_VI_INF);  CHKERRQ(ierr);
  if (sys.nonlinear_solver->bounds != NULL)
    sys.nonlinear_solver->bounds(XL,XU,sys);
  else if (sys.nonlinear_solver->bounds_object != NULL)
    sys.nonlinear_solver->bounds_object->bounds(XL,XU, sys);
  else
    SETERRQ(((PetscObject)dm)->comm, PETSC_ERR_ARG_WRONG, "No bounds calculation in this libMesh object");

  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "DMCreateGlobalVector_libMesh"
static PetscErrorCode DMCreateGlobalVector_libMesh(DM dm, Vec *x)
{
  PetscFunctionBegin;
  PetscErrorCode ierr;
  DM_libMesh     *dlm = (DM_libMesh *)(dm->data);
  PetscBool eq;

  ierr = PetscObjectTypeCompare((PetscObject)dm, DMLIBMESH, &eq); CHKERRQ(ierr);

  if (!eq)
    SETERRQ2(((PetscObject)dm)->comm, PETSC_ERR_ARG_WRONG, "DM of type %s, not of type %s", ((PetscObject)dm)->type, DMLIBMESH);

  if (!dlm->sys)
    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONGSTATE, "No libMesh system set for DM_libMesh");

  NumericVector<Number>* nv = (dlm->sys->solution).get();
  PetscVector<Number>*   pv = dynamic_cast<PetscVector<Number>*>(nv);
  Vec                    v  = pv->vec();
  /* Unfortunately, currently this does not produce a ghosted vector, so nonlinear subproblem solves aren't going to be easily available.
     Should work fine for getting vectors out for linear subproblem solvers. */
  if(dlm->embedding) {
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
  ierr = PetscObjectCompose((PetscObject)*x,"DM",(PetscObject)dm); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}




#undef __FUNCT__
#define __FUNCT__ "DMCreateMatrix_libMesh"
#if PETSC_VERSION_LT(3,5,0)
static PetscErrorCode DMCreateMatrix_libMesh(DM dm, const MatType, Mat *A)
#else
static PetscErrorCode DMCreateMatrix_libMesh(DM dm, Mat *A)
#endif
{
  PetscFunctionBegin;
  PetscErrorCode ierr;
  DM_libMesh     *dlm = (DM_libMesh *)(dm->data);
  PetscBool eq;

  ierr = PetscObjectTypeCompare((PetscObject)dm, DMLIBMESH, &eq); CHKERRQ(ierr);

  if (!eq)
    SETERRQ2(((PetscObject)dm)->comm, PETSC_ERR_ARG_WRONG, "DM of type %s, not of type %s", ((PetscObject)dm)->type, DMLIBMESH);

  if (!dlm->sys)
    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONGSTATE, "No libMesh system set for DM_libMesh");

  *A = (dynamic_cast<PetscMatrix<Number>*>(dlm->sys->matrix))->mat();
  ierr = PetscObjectReference((PetscObject)(*A)); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "DMView_libMesh"
static PetscErrorCode  DMView_libMesh(DM dm, PetscViewer viewer)
{
  PetscErrorCode ierr;
  PetscBool isascii;
  const char *name, *prefix;
  DM_libMesh *dlm = (DM_libMesh*)dm->data;
  PetscFunctionBegin;
  ierr = PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii); CHKERRQ(ierr);
  if(isascii) {
    ierr = PetscObjectGetName((PetscObject)dm, &name);     CHKERRQ(ierr);
    ierr = PetscObjectGetOptionsPrefix((PetscObject)dm, &prefix); CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer, "DM libMesh with name %s and prefix %s\n", name, prefix); CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer, "blocks:", name, prefix); CHKERRQ(ierr);
    std::map<std::string,unsigned int>::iterator bit = dlm->blockids->begin();
    std::map<std::string,unsigned int>::const_iterator bend = dlm->blockids->end();
    for(; bit != bend; ++bit) {
      ierr = PetscViewerASCIIPrintf(viewer, "(%s,%D) ", bit->first.c_str(), bit->second); CHKERRQ(ierr);
    }
    ierr = PetscViewerASCIIPrintf(viewer, "\n"); CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer, "variables:", name, prefix); CHKERRQ(ierr);
    std::map<std::string,unsigned int>::iterator vit = dlm->varids->begin();
    std::map<std::string,unsigned int>::const_iterator vend = dlm->varids->end();
    for(; vit != vend; ++vit) {
      ierr = PetscViewerASCIIPrintf(viewer, "(%s,%D) ", vit->first.c_str(), vit->second); CHKERRQ(ierr);
    }
    ierr = PetscViewerASCIIPrintf(viewer, "\n"); CHKERRQ(ierr);
    if(dlm->decomposition_type == DMLIBMESH_NO_DECOMPOSITION) {
      ierr = PetscViewerASCIIPrintf(viewer, "No decomposition\n"); CHKERRQ(ierr);
    }
    else {
      if(dlm->decomposition_type == DMLIBMESH_FIELD_DECOMPOSITION) {
	ierr = PetscViewerASCIIPrintf(viewer, "Field decomposition by variable: "); CHKERRQ(ierr);
      }
      else if(dlm->decomposition_type == DMLIBMESH_DOMAIN_DECOMPOSITION) {
	ierr = PetscViewerASCIIPrintf(viewer, "Domain decomposition by block: "); CHKERRQ(ierr);
      }
      else SETERRQ1(((PetscObject)dm)->comm, PETSC_ERR_PLIB, "Unexpected decomposition type: %D", dlm->decomposition_type);
      /* FIX: decompositions might have different sizes and components on different ranks. */
      for(unsigned int d = 0; d < dlm->decomposition->size(); ++d) {
	std::set<unsigned int>::iterator dbegin  = (*dlm->decomposition)[d].begin();
	std::set<unsigned int>::iterator dit     = (*dlm->decomposition)[d].begin();
	std::set<unsigned int>::iterator dend    = (*dlm->decomposition)[d].end();
	for(; dit != dend; ++dit) {
	  if(dit != dbegin) {
	    ierr = PetscViewerASCIIPrintf(viewer, ","); CHKERRQ(ierr);
	  }
	  ierr = PetscViewerASCIIPrintf(viewer, "%D", *dit); CHKERRQ(ierr);
	}
	ierr = PetscViewerASCIIPrintf(viewer, ";"); CHKERRQ(ierr);
      }
      ierr = PetscViewerASCIIPrintf(viewer, "\n"); CHKERRQ(ierr);
    }
  }

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMSetUp_libMesh"
static PetscErrorCode  DMSetUp_libMesh(DM dm)
{
  PetscFunctionBegin;
  PetscErrorCode ierr;
  DM_libMesh     *dlm = (DM_libMesh *)(dm->data);
  PetscBool eq;

  ierr = PetscObjectTypeCompare((PetscObject)dm, DMLIBMESH, &eq); CHKERRQ(ierr);

  if (!eq)
    SETERRQ2(((PetscObject)dm)->comm, PETSC_ERR_ARG_WRONG, "DM of type %s, not of type %s", ((PetscObject)dm)->type, DMLIBMESH);

  if (!dlm->sys)
    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONGSTATE, "No libMesh system set for DM_libMesh");
  /*
     Do not evaluate function, Jacobian or bounds for an embedded DM -- the subproblem might not have enough information for that.
  */
  if(!dlm->embedding) {
#if PETSC_RELEASE_LESS_THAN(3,3,1)
    ierr = DMSetFunction(dm, DMlibMeshFunction); CHKERRQ(ierr);
    ierr = DMSetJacobian(dm, DMlibMeshJacobian); CHKERRQ(ierr);
#else
    ierr = DMSNESSetFunction(dm, SNESFunction_DMlibMesh, (void*)dm); CHKERRQ(ierr);
    ierr = DMSNESSetJacobian(dm, SNESJacobian_DMlibMesh, (void*)dm); CHKERRQ(ierr);
#endif
    if (dlm->sys->nonlinear_solver->bounds || dlm->sys->nonlinear_solver->bounds_object)
      ierr = DMSetVariableBounds(dm, DMVariableBounds_libMesh); CHKERRQ(ierr);
  }
  else {
    /*
       Fow now we don't implement even these, although a linear "Dirichlet" subproblem is well-defined.
       Creating the submatrix, however, might require extracting the submatrix preallocation from an unassembled matrix.
    */
      dm->ops->createglobalvector = 0;
      dm->ops->creatematrix = 0;
  }
  PetscFunctionReturn(0);
}




#undef __FUNCT__
#define __FUNCT__ "DMDestroy_libMesh"
static PetscErrorCode  DMDestroy_libMesh(DM dm)
{
  DM_libMesh *dlm = (DM_libMesh*)(dm->data);
  PetscErrorCode ierr;
  PetscFunctionBegin;
  delete dlm->varids;
  delete dlm->varnames;
  delete dlm->blockids;
  delete dlm->blocknames;
  delete dlm->decomposition;
  ierr = ISDestroy(&dlm->embedding); CHKERRQ(ierr);
  ierr = PetscFree(dm->data); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "DMCreateLibMesh"
PetscErrorCode  DMCreateLibMesh(MPI_Comm comm, NonlinearImplicitSystem& sys, DM *dm)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = DMCreate(comm, dm);           CHKERRQ(ierr);
  ierr = DMSetType(*dm, DMLIBMESH);    CHKERRQ(ierr);
  ierr = DMLibMeshSetSystem(*dm, sys); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

EXTERN_C_BEGIN
#undef __FUNCT__
#define __FUNCT__ "DMCreate_libMesh"
PetscErrorCode  DMCreate_libMesh(DM dm)
{
  PetscErrorCode ierr;
  DM_libMesh     *dlm;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
#if PETSC_RELEASE_LESS_THAN(3,5,0)
  ierr = PetscNewLog(dm,DM_libMesh,&dlm);CHKERRQ(ierr);
#else
  ierr = PetscNewLog(dm,&dlm);CHKERRQ(ierr);
#endif
  dm->data = dlm;

  dlm->varids     = new(std::map<std::string, unsigned int>);
  dlm->blockids   = new(std::map<std::string, unsigned int>);
  dlm->varnames   = new(std::map<unsigned int, std::string>);
  dlm->blocknames = new(std::map<unsigned int, std::string>);
  dlm->decomposition   = PETSC_NULL;
  dlm->decomposition_type  = DMLIBMESH_NO_DECOMPOSITION;

  dm->ops->createglobalvector = DMCreateGlobalVector_libMesh;
  dm->ops->createlocalvector  = 0; // DMCreateLocalVector_libMesh;
  dm->ops->getcoloring        = 0; // DMGetColoring_libMesh;
  dm->ops->creatematrix       = DMCreateMatrix_libMesh;
  dm->ops->createinterpolation= 0; // DMCreateInterpolation_libMesh;

  dm->ops->refine             = 0; // DMRefine_libMesh;
  dm->ops->coarsen            = 0; // DMCoarsen_libMesh;
  dm->ops->getinjection       = 0; // DMGetInjection_libMesh;
  dm->ops->getaggregates      = 0; // DMGetAggregates_libMesh;

#if PETSC_RELEASE_LESS_THAN(3,3,1)
  dm->ops->createfielddecompositiondm  = DMCreateFieldDecompositionDM_libMesh;
  dm->ops->createdomaindecompositiondm = DMCreateDomainDecompositionDM_libMesh;
#endif
  dm->ops->createfielddecomposition    = DMCreateFieldDecomposition_libMesh;
  dm->ops->createdomaindecomposition   = DMCreateDomainDecomposition_libMesh;

  dm->ops->destroy            = DMDestroy_libMesh;
  dm->ops->view               = DMView_libMesh;
  dm->ops->setfromoptions     = 0; // DMSetFromOptions_libMesh;
  dm->ops->setup              = DMSetUp_libMesh;

  PetscFunctionReturn(0);
}
EXTERN_C_END


#endif // #if !PETSC_VERSION_LESS_THAN(3,3,0)
