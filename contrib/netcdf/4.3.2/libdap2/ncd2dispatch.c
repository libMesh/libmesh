/*********************************************************************
 *   Copyright 1993, UCAR/Unidata
 *   See netcdf/COPYRIGHT file for copying and redistribution conditions.
 *********************************************************************/

#include "ncdap.h"
#include "ncd2dispatch.h"

#ifdef HAVE_GETRLIMIT
#  ifdef HAVE_SYS_RESOURCE_H
#    include <sys/time.h>
#  endif
#  ifdef HAVE_SYS_RESOURCE_H
#    include <sys/resource.h>
#  endif
#endif

#define getncid(drno) (((NC*)drno)->ext_ncid)

/* Define the set of protocols known to be constrainable */
static char* constrainableprotocols[] = {"http", "https",NULL};

static int ncd2initialized = 0;

size_t dap_one[NC_MAX_VAR_DIMS];
size_t dap_zero[NC_MAX_VAR_DIMS];

/* Forward */
static NCerror buildncstructures(NCDAPCOMMON*);
static NCerror builddims(NCDAPCOMMON*);
static char* getdefinename(CDFnode* node);
static NCerror buildvars(NCDAPCOMMON*);
static NCerror buildglobalattrs(NCDAPCOMMON*, CDFnode* root);
static NCerror buildattribute(NCDAPCOMMON*, NCattribute*, nc_type, int);
static void computedimindexanon(CDFnode* dim, CDFnode* var);
static void replacedims(NClist* dims);
static int equivalentdim(CDFnode* basedim, CDFnode* dupdim);
static NCerror addstringdims(NCDAPCOMMON*);
static NCerror defrecorddim(NCDAPCOMMON*);
static NCerror defseqdims(NCDAPCOMMON*);
static NCerror showprojection(NCDAPCOMMON*, CDFnode* var);
static NCerror getseqdimsize(NCDAPCOMMON*, CDFnode* seq, size_t* sizep);
static NCerror makeseqdim(NCDAPCOMMON*, CDFnode* seq, size_t count, CDFnode** sqdimp);
static NCerror countsequence(NCDAPCOMMON*, CDFnode* xseq, size_t* sizep);
static NCerror freeNCDAPCOMMON(NCDAPCOMMON*);
static NCerror fetchtemplatemetadata(NCDAPCOMMON*);
static int fieldindex(CDFnode* parent, CDFnode* child);
static NCerror computeseqcountconstraints(NCDAPCOMMON*, CDFnode*, NCbytes*);
static void computeseqcountconstraintsr(NCDAPCOMMON*, CDFnode*, CDFnode**);
static void estimatevarsizes(NCDAPCOMMON*);
static NCerror fetchconstrainedmetadata(NCDAPCOMMON*);
static NCerror suppressunusablevars(NCDAPCOMMON*);
static NCerror fixzerodims(NCDAPCOMMON*);
static void applyclientparamcontrols(NCDAPCOMMON*);
static NCerror applyclientparams(NCDAPCOMMON*);

/**************************************************/

static int
NCD2_create(const char *path, int cmode,
           size_t initialsz, int basepe, size_t *chunksizehintp,
	   int use_parallel, void* mpidata,
           NC_Dispatch*,NC* ncp);

static int NCD2_redef(int ncid);
static int NCD2__enddef(int ncid, size_t h_minfree, size_t v_align, size_t v_minfree, size_t r_align);
static int NCD2_sync(int ncid);
static int NCD2_abort(int ncid);

static int NCD2_put_vara(int ncid, int varid,
	    const size_t *start, const size_t *edges0,
            const void *value0,
	    nc_type memtype);

static int NCD2_get_vara(int ncid, int varid,
	    const size_t *start, const size_t *edges,
            void *value,
	    nc_type memtype);

static int NCD2_put_vars(int ncid, int varid,
	    const size_t *start, const size_t *edges, const ptrdiff_t* stride,
            const void *value0, nc_type memtype);

static int NCD2_get_vars(int ncid, int varid,
	    const size_t *start, const size_t *edges, const ptrdiff_t* stride,
            void *value, nc_type memtype);

static NC_Dispatch NCD2_dispatch_base = {

NC_DISPATCH_NC3 | NC_DISPATCH_NCD,

NCD2_create,
NCD2_open,

NCD2_redef,
NCD2__enddef,
NCD2_sync,
NCD2_abort,
NCD2_close,
NULL, /*set_fill*/
NULL, /*inq_base_pe*/
NULL, /*set_base_pe*/
NULL, /*inq_format*/
NCD2_inq_format_extended, /*inq_format_extended*/

NULL, /*inq*/
NULL, /*inq_type*/

NULL, /*def_dim*/
NULL, /*inq_dimid*/
NULL, /*inq_dim*/
NULL, /*inq_unlimdim*/
NULL, /*rename_dim*/

NULL, /*inq_att*/
NULL, /*inq_attid*/
NULL, /*inq_attname*/
NULL, /*rename_att*/
NULL, /*del_att*/
NULL, /*get_att*/
NULL, /*put_att*/

NULL, /*def_var*/
NULL, /*inq_varid*/
NULL, /*rename_var*/
NCD2_get_vara,
NCD2_put_vara,
NCD2_get_vars,
NCD2_put_vars,
NCDEFAULT_get_varm,
NCDEFAULT_put_varm,

NULL, /*inq_var_all*/

NULL, /*var_par_access*/

#ifdef USE_NETCDF4
NULL, /*show_metadata*/
NULL, /*inq_unlimdims*/
NULL, /*inq_ncid*/
NULL, /*inq_grps*/
NULL, /*inq_grpname*/
NULL, /*inq_grpname_full*/
NULL, /*inq_grp_parent*/
NULL, /*inq_grp_full_ncid*/
NULL, /*inq_varids*/
NULL, /*inq_dimids*/
NULL, /*inq_typeids*/
NULL, /*inq_type_equal*/
NULL, /*def_grp*/
NULL, /*rename_grp*/
NULL, /*inq_user_type*/
NULL, /*inq_typeid*/

NULL, /*def_compound*/
NULL, /*insert_compound*/
NULL, /*insert_array_compound*/
NULL, /*inq_compound_field*/
NULL, /*inq_compound_fieldindex*/
NULL, /*def_vlen*/
NULL, /*put_vlen_element*/
NULL, /*get_vlen_element*/
NULL, /*def_enum*/
NULL, /*insert_enum*/
NULL, /*inq_enum_member*/
NULL, /*inq_enum_ident*/
NULL, /*def_opaque*/
NULL, /*def_var_deflate*/
NULL, /*def_var_fletcher32*/
NULL, /*def_var_chunking*/
NULL, /*def_var_fill*/
NULL, /*def_var_endian*/
NULL, /*set_var_chunk_cache*/
NULL, /*get_var_chunk_cache*/

#endif /*USE_NETCDF4*/

};

NC_Dispatch* NCD2_dispatch_table = NULL; /* moved here from ddispatch.c */

static NC_Dispatch NCD2_dispatcher; /* overlay result */

int
NCD2_initialize(void)
{
    int i;
    /* Create our dispatch table as the merge of NCD2 table and NCSUBSTRATE */
    /* watch the order because we want NCD2 to overwrite NCSUBSTRATE */
    NC_dispatch_overlay(&NCD2_dispatch_base, NCSUBSTRATE_dispatch_table, &NCD2_dispatcher);    
    NCD2_dispatch_table = &NCD2_dispatcher;
    /* Local Initialization */
    compute_nccalignments();
    for(i=0;i<NC_MAX_VAR_DIMS;i++) {
	dap_one[i] = 1;
	dap_zero[i] = 0;
    }
    ncd2initialized = 1;
#ifdef DEBUG
    /* force logging to go to stderr */
    nclogclose();
    if(nclogopen(NULL))
        ncsetlogging(1); /* turn it on */
#endif
    return NC_NOERR;
}

static int
NCD2_redef(int ncid)
{
    return (NC_EPERM);
}

static int
NCD2__enddef(int ncid, size_t h_minfree, size_t v_align, size_t v_minfree, size_t r_align)
{
    return (NC_EPERM);
}

static int
NCD2_sync(int ncid)
{
    return (NC_EINVAL);
}

static int
NCD2_abort(int ncid)
{
    return NCD2_close(ncid);
}

static int
NCD2_create(const char *path, int cmode,
           size_t initialsz, int basepe, size_t *chunksizehintp,
	   int use_parallel, void* mpidata,
           NC_Dispatch* dispatch, NC* ncp)
{
   return NC_EPERM;
}

static int
NCD2_put_vara(int ncid, int varid,
	    const size_t *start, const size_t *edges,
            const void *value,
	    nc_type memtype)
{
    return NC_EPERM;
}

static int
NCD2_get_vara(int ncid, int varid,
	    const size_t *start, const size_t *edges,
            void *value,
	    nc_type memtype)
{
    int stat = nc3d_getvarx(ncid, varid, start, edges, nc_ptrdiffvector1, value,memtype);
    return stat;
}

static int
NCD2_put_vars(int ncid, int varid,
	    const size_t *start, const size_t *edges, const ptrdiff_t* stride,
            const void *value0, nc_type memtype)
{
    return NC_EPERM;
}

static int
NCD2_get_vars(int ncid, int varid,
	    const size_t *start, const size_t *edges, const ptrdiff_t* stride,
            void *value, nc_type memtype)
{
    int stat = nc3d_getvarx(ncid, varid, start, edges, stride, value, memtype);
    return stat;
}

/* See ncd2dispatch.c for other version */
int
NCD2_open(const char * path, int mode,
               int basepe, size_t *chunksizehintp,
 	       int useparallel, void* mpidata,
               NC_Dispatch* dispatch, NC* drno)
{
    NCerror ncstat = NC_NOERR;
    OCerror ocstat = OC_NOERR;
    NCDAPCOMMON* dapcomm = NULL;
    const char* value;

    if(path == NULL)
	return NC_EDAPURL;
    if(dispatch == NULL) PANIC("NC3D_open: no dispatch table");

    /* Setup our NC and NCDAPCOMMON state*/

    dapcomm = (NCDAPCOMMON*)calloc(1,sizeof(NCDAPCOMMON));
    if(dapcomm == NULL) {ncstat = NC_ENOMEM; goto done;}

    NCD2_DATA_SET(drno,dapcomm);
    drno->int_ncid = nc__pseudofd(); /* create a unique id */
    dapcomm->controller = (NC*)drno;

    dapcomm->cdf.separator = ".";
    dapcomm->cdf.smallsizelimit = DFALTSMALLLIMIT;
    dapcomm->cdf.cache = createnccache();

#ifdef HAVE_GETRLIMIT
    { struct rlimit rl;
      if(getrlimit(RLIMIT_NOFILE, &rl) >= 0) {
	dapcomm->cdf.cache->cachecount = (size_t)(rl.rlim_cur / 2);
      }
    }
#endif

#ifdef OCCOMPILEBYDEFAULT
    /* set the compile flag by default */
    dapcomm->oc.rawurltext = (char*)emalloc(strlen(path)+strlen("[compile]")+1);
    strcpy(dapcomm->oc.rawurltext,"[compile]");
    strcat(dapcomm->oc.rawurltext, path);    
#else
    dapcomm->oc.rawurltext = strdup(path);
#endif

    ncuriparse(dapcomm->oc.rawurltext,&dapcomm->oc.url);

    /* parse the client parameters */
    ncuridecodeparams(dapcomm->oc.url);

    if(!constrainable(dapcomm->oc.url))
	SETFLAG(dapcomm->controls,NCF_UNCONSTRAINABLE);

#ifdef COLUMBIA_HACK
    {
	const char* p;
	/* Does this url look like it is from columbia? */
	if(dapcomm->oc.url->host != NULL) {
	    for(p=dapcomm->oc.url->host;*p;p++) {
	        if(strncmp(p,COLUMBIA_HACK,strlen(COLUMBIA_HACK))==0)
		    SETFLAG(dapcomm->controls,NCF_COLUMBIA);
	    }
	}
    } 
#endif

    /* fail if we are unconstrainable but have constraints */
    if(FLAGSET(dapcomm->controls,NCF_UNCONSTRAINABLE)) {
	if(dapcomm->oc.url->constraint != NULL) {
	    nclog(NCLOGWARN,"Attempt to constrain an unconstrainable data source: %s",
		   dapcomm->oc.url->constraint);
	    ncstat = THROW(NC_EDAPCONSTRAINT);
	    goto done;
	}
    }

    /* Use libsrc code for storing metadata */
    {
	char tmpname[32];

        /* Create fake file name: exact name must be unique,
           but is otherwise irrelevant because we are using NC_DISKLESS
        */
        snprintf(tmpname,sizeof(tmpname),"%d",drno->int_ncid);

        /* Now, use the file to create the netcdf file; force classic.  */
        ncstat = nc_create(tmpname,NC_DISKLESS|NC_CLASSIC_MODEL,&drno->substrate);
        if(ncstat != NC_NOERR) {THROWCHK(ncstat); goto done;}
    }

    /* Avoid fill */
    nc_set_fill(drno->substrate,NC_NOFILL,NULL);

    dapcomm->oc.dapconstraint = (DCEconstraint*)dcecreate(CES_CONSTRAINT);
    dapcomm->oc.dapconstraint->projections = nclistnew();
    dapcomm->oc.dapconstraint->selections = nclistnew();
    
     /* Parse constraints to make sure they are syntactically correct */
     ncstat = parsedapconstraints(dapcomm,dapcomm->oc.url->constraint,dapcomm->oc.dapconstraint);
     if(ncstat != NC_NOERR) {THROWCHK(ncstat); goto done;}

    /* Construct a url for oc minus any constraint and params*/
    dapcomm->oc.urltext = ncuribuild(dapcomm->oc.url,NULL,NULL,
				      (NCURISTD ^ NCURICONSTRAINTS));

    /* Pass to OC */
    ocstat = oc_open(dapcomm->oc.urltext,&dapcomm->oc.conn);
    if(ocstat != OC_NOERR) {THROWCHK(ocstat); goto done;}

#ifdef DEBUG1
    (void)oc_trace_curl(dapcomm->oc.conn);
#endif

    nullfree(dapcomm->oc.urltext); /* clean up */
    dapcomm->oc.urltext = NULL;

    /* process control client parameters */
    applyclientparamcontrols(dapcomm);

    /* Turn on logging; only do this after oc_open*/
    if((value = dapparamvalue(dapcomm,"log")) != NULL) {
	ncloginit();
        if(nclogopen(value))
	    ncsetlogging(1);
	ocloginit();
        if(oclogopen(value))
	    ocsetlogging(1);
    }

    /* fetch and build the unconstrained DDS for use as
       template */
    ncstat = fetchtemplatemetadata(dapcomm);
    if(ncstat != NC_NOERR) goto done;

    /* Operations on the template tree */

    /* Accumulate useful nodes sets  */
    ncstat = computecdfnodesets(dapcomm,dapcomm->cdf.fullddsroot->tree);
    if(ncstat) {THROWCHK(ncstat); goto done;}

    /* Define the dimsettrans list */
    ncstat = definedimsettrans(dapcomm,dapcomm->cdf.fullddsroot->tree);
    if(ncstat) {THROWCHK(ncstat); goto done;}

    /* Mark the nodes of the template that are eligible for prefetch */
    ncstat = markprefetch(dapcomm);

    /* fetch and build the constrained DDS */
    ncstat = fetchconstrainedmetadata(dapcomm);
    if(ncstat != NC_NOERR) goto done;

#ifdef DEBUG2
fprintf(stderr,"constrained dds: %s\n",dumptree(dapcomm->cdf.ddsroot));
#endif

    /* Operations on the constrained tree */

    /* Accumulate useful nodes sets  */
    ncstat = computecdfnodesets(dapcomm,dapcomm->cdf.ddsroot->tree);
    if(ncstat) {THROWCHK(ncstat); goto done;}

    /* Fix grids */
    ncstat = fixgrids(dapcomm);
    if(ncstat) {THROWCHK(ncstat); goto done;}

    /* Locate and mark usable sequences */
    ncstat = sequencecheck(dapcomm);
    if(ncstat) {THROWCHK(ncstat); goto done;}

    /* suppress variables not in usable sequences */
    ncstat = suppressunusablevars(dapcomm);
    if(ncstat) {THROWCHK(ncstat); goto done;}

    /* apply client parameters */
    ncstat = applyclientparams(dapcomm);
    if(ncstat) {THROWCHK(ncstat); goto done;}

    /* Add (as needed) string dimensions*/
    ncstat = addstringdims(dapcomm);
    if(ncstat) {THROWCHK(ncstat); goto done;}

    if(nclistlength(dapcomm->cdf.ddsroot->tree->seqnodes) > 0) {
	/* Build the sequence related dimensions */
        ncstat = defseqdims(dapcomm);
        if(ncstat) {THROWCHK(ncstat); goto done;}
    }

    /* Define the dimsetplus and dimsetall lists */
    ncstat = definedimsets(dapcomm,dapcomm->cdf.ddsroot->tree);
    if(ncstat) {THROWCHK(ncstat); goto done;}

    /* Re-compute the dimension names*/
    ncstat = computecdfdimnames(dapcomm);
    if(ncstat) {THROWCHK(ncstat); goto done;}

    /* Deal with zero size dimensions */
    ncstat = fixzerodims(dapcomm);
    if(ncstat) {THROWCHK(ncstat); goto done;}

    /* Attempt to use the DODS_EXTRA info to turn
       one of the dimensions into unlimited.
       Assume computecdfdimnames34 has already been called.
    */
    ncstat = defrecorddim(dapcomm);
    if(ncstat) {THROWCHK(ncstat); goto done;}
    if(dapcomm->cdf.recorddimname != NULL
       && nclistlength(dapcomm->cdf.ddsroot->tree->seqnodes) > 0) {
	/*nclog(NCLOGWARN,"unlimited dimension specified, but sequences exist in DDS");*/
	PANIC("unlimited dimension specified, but sequences exist in DDS");	
    }

    /* Re-compute the var names*/
    ncstat = computecdfvarnames(dapcomm,
				 dapcomm->cdf.ddsroot,
				 dapcomm->cdf.ddsroot->tree->varnodes);
    if(ncstat) {THROWCHK(ncstat); goto done;}

    /* Transfer data from the unconstrained DDS data to the unconstrained DDS */
    ncstat = dimimprint(dapcomm);
    if(ncstat) goto done;

    /* Process the constraints to map to the constrained CDF tree */
    /* (must follow fixgrids3 */
    ncstat = mapconstraints(dapcomm->oc.dapconstraint,dapcomm->cdf.ddsroot);
    if(ncstat != NC_NOERR) goto done;

    /* Canonicalize the constraint */
    ncstat = fixprojections(dapcomm->oc.dapconstraint->projections);
    if(ncstat != NC_NOERR) {THROWCHK(ncstat); goto done;}

    /* Fill in segment information */
    ncstat = qualifyconstraints(dapcomm->oc.dapconstraint);
    if(ncstat != NC_NOERR) goto done;

    /* Accumulate set of variables in the constraint's projections */
    ncstat = computeprojectedvars(dapcomm,dapcomm->oc.dapconstraint);
    if(ncstat) {THROWCHK(ncstat); goto done;}

    /* using the modified constraint, rebuild the constraint string */
    if(FLAGSET(dapcomm->controls,NCF_UNCONSTRAINABLE)) {
	/* ignore all constraints */
	dapcomm->oc.urltext = ncuribuild(dapcomm->oc.url,NULL,NULL,0);
    } else {
	char* constraintstring = buildconstraintstring(dapcomm->oc.dapconstraint);
        ncurisetconstraints(dapcomm->oc.url,constraintstring);
	nullfree(constraintstring);
        dapcomm->oc.urltext = ncuribuild(dapcomm->oc.url,NULL,NULL,NCURICONSTRAINTS);
    }

#ifdef DEBUG
fprintf(stderr,"ncdap3: final constraint: %s\n",dapcomm->oc.url->constraint);
#endif

    /* Estimate the variable sizes */
    estimatevarsizes(dapcomm);

    /* Build the meta data */
    ncstat = buildncstructures(dapcomm);
    if(ncstat != NC_NOERR) {THROWCHK(ncstat); goto done;}

    /* Explicitly do not call enddef because it will complain
       about variables that are too large.
    */
#if 0
    ncstat = nc_endef(drno->substrate,NC_NOFILL,NULL);
    if(ncstat != NC_NOERR && ncstat != NC_EVARSIZE)
        {THROWCHK(ncstat); goto done;}
#endif
	
    {
        NC* ncsub;
        NC* drno = dapcomm->controller;
	CDFnode* unlimited = dapcomm->cdf.recorddim;
        /* (for now) break abstractions*/
	NC3_INFO* nc3i;

        /* get the id for the substrate */
        ncstat = NC_check_id(drno->substrate,&ncsub);
        if(ncstat != NC_NOERR) {THROWCHK(ncstat); goto done;}
	nc3i = (NC3_INFO*)ncsub->dispatchdata;

        if(unlimited != NULL) {
            /* Set the effective size of UNLIMITED */
            NC_set_numrecs(nc3i,unlimited->dim.declsize);
        }

        /* Pretend the substrate is read-only */
	NC_set_readonly(nc3i);
	
    }

    /* Do any necessary data prefetch */
    if(FLAGSET(dapcomm->controls,NCF_PREFETCH)
       && FLAGSET(dapcomm->controls,NCF_PREFETCH_EAGER)) {
        ncstat = prefetchdata(dapcomm);
        if(ncstat != NC_NOERR) {
            del_from_NCList((NC*)drno); /* undefine here */
	    {THROWCHK(ncstat); goto done;}
	}
    }

    return ncstat;

done:
    if(drno != NULL) NCD2_close(drno->ext_ncid);
    if(ocstat != OC_NOERR) ncstat = ocerrtoncerr(ocstat);
    return THROW(ncstat);
}

int
NCD2_close(int ncid)
{
    NC* drno;
    NCDAPCOMMON* dapcomm;
    int ncstatus = NC_NOERR;

    ncstatus = NC_check_id(ncid, (NC**)&drno); 
    if(ncstatus != NC_NOERR) return THROW(ncstatus);
    dapcomm = (NCDAPCOMMON*)drno->dispatchdata;

    /* We call abort rather than close to avoid
       trying to write anything or try to pad file length
     */
    ncstatus = nc_abort(drno->substrate);

    /* clean NC* */
    freeNCDAPCOMMON(dapcomm);

    return THROW(ncstatus);
}

/**************************************************/

static NCerror
buildncstructures(NCDAPCOMMON* dapcomm)
{
    NCerror ncstat = NC_NOERR;
    CDFnode* dds = dapcomm->cdf.ddsroot;
    NC* ncsub;

    ncstat = NC_check_id(dapcomm->controller->substrate,&ncsub);
    if(ncstat != NC_NOERR) goto done;

    ncstat = buildglobalattrs(dapcomm,dds);
    if(ncstat != NC_NOERR) goto done;

    ncstat = builddims(dapcomm);
    if(ncstat != NC_NOERR) goto done;

    ncstat = buildvars(dapcomm);
    if(ncstat != NC_NOERR) goto done;

done:
    return THROW(ncstat);
}

static NCerror
builddims(NCDAPCOMMON* dapcomm)
{
    int i;
    NCerror ncstat = NC_NOERR;
    int dimid;
    NClist* dimset = NULL;
    NC* drno = dapcomm->controller;
    NC* ncsub;
    char* definename;

    /* collect all dimensions from variables */
    dimset = dapcomm->cdf.ddsroot->tree->dimnodes;

    /* Sort by fullname just for the fun of it */
    for(;;) {
	int last = nclistlength(dimset) - 1;
	int swap = 0;
        for(i=0;i<last;i++) {
	    CDFnode* dim1 = (CDFnode*)nclistget(dimset,i);
	    CDFnode* dim2 = (CDFnode*)nclistget(dimset,i+1);
   	    if(strcmp(dim1->ncfullname,dim2->ncfullname) > 0) {
		nclistset(dimset,i,(void*)dim2);
		nclistset(dimset,i+1,(void*)dim1);
		swap = 1;
		break;
	    }
	}
	if(!swap) break;
    }

    /* Define unlimited only if needed */ 
    if(dapcomm->cdf.recorddim != NULL) {
	CDFnode* unlimited = dapcomm->cdf.recorddim;
	definename = getdefinename(unlimited);
        ncstat = nc_def_dim(drno->substrate,
			definename,
			NC_UNLIMITED,
			&unlimited->ncid);
	nullfree(definename);
        if(ncstat != NC_NOERR) {THROWCHK(ncstat); goto done;}

        /* get the id for the substrate */
        ncstat = NC_check_id(drno->substrate,&ncsub);
        if(ncstat != NC_NOERR) {THROWCHK(ncstat); goto done;}
#if 0
	nc3sub = (NC3_INFO*)&ncsub->dispatchdata;
        /* Set the effective size of UNLIMITED;
           note that this cannot easily be done thru the normal API.*/
        NC_set_numrecs(nc3sub,unlimited->dim.declsize);
#endif

    }

    for(i=0;i<nclistlength(dimset);i++) {
	CDFnode* dim = (CDFnode*)nclistget(dimset,i);
        if(dim->dim.basedim != NULL) continue; /* handle below */
	if(DIMFLAG(dim,CDFDIMRECORD)) continue; /* defined above */
#ifdef DEBUG1
fprintf(stderr,"define: dim: %s=%ld\n",dim->ncfullname,(long)dim->dim.declsize);
#endif
	definename = getdefinename(dim);
        ncstat = nc_def_dim(drno->substrate,definename,dim->dim.declsize,&dimid);
        if(ncstat != NC_NOERR) {
	    THROWCHK(ncstat); goto done;
	}
	nullfree(definename);
        dim->ncid = dimid;
    }

    /* Make all duplicate dims have same dimid as basedim*/
    /* (see computecdfdimnames)*/
    for(i=0;i<nclistlength(dimset);i++) {
	CDFnode* dim = (CDFnode*)nclistget(dimset,i);
        if(dim->dim.basedim != NULL) {
	    dim->ncid = dim->dim.basedim->ncid;
	}
    }
done:
    nclistfree(dimset);
    return THROW(ncstat);
}

/* Simultaneously build any associated attributes*/
/* and any necessary pseudo-dimensions for string types*/
static NCerror
buildvars(NCDAPCOMMON* dapcomm)
{
    int i,j;
    NCerror ncstat = NC_NOERR;
    int varid;
    NClist* varnodes = dapcomm->cdf.ddsroot->tree->varnodes;
    NC* drno = dapcomm->controller;
    char* definename;

    ASSERT((varnodes != NULL));
    for(i=0;i<nclistlength(varnodes);i++) {
	CDFnode* var = (CDFnode*)nclistget(varnodes,i);
        int dimids[NC_MAX_VAR_DIMS];
	unsigned int ncrank;
        NClist* vardims = NULL;

	if(var->invisible) continue;
	if(var->array.basevar != NULL) continue;

#ifdef DEBUG1
fprintf(stderr,"buildvars.candidate=|%s|\n",var->ncfullname);
#endif

	vardims = var->array.dimsetall;
	ncrank = nclistlength(vardims);
	if(ncrank > 0) {
            for(j=0;j<ncrank;j++) {
                CDFnode* dim = (CDFnode*)nclistget(vardims,j);
                dimids[j] = dim->ncid;
 	    }
        }   



	definename = getdefinename(var);

#ifdef DEBUG1
fprintf(stderr,"define: var: %s/%s",
		definename,var->ocname);
if(ncrank > 0) {
int k;
for(k=0;k<ncrank;k++) {
CDFnode* dim = (CDFnode*)nclistget(vardims,k);
fprintf(stderr,"[%ld]",dim->dim.declsize);
 }
 }
fprintf(stderr,"\n");
#endif
        ncstat = nc_def_var(drno->substrate,
		        definename,
                        var->externaltype,
                        ncrank,
                        (ncrank==0?NULL:dimids),
                        &varid);
	nullfree(definename);
        if(ncstat != NC_NOERR) {
	    THROWCHK(ncstat);
	    goto done;
	}
        var->ncid = varid;
	if(var->attributes != NULL) {
	    for(j=0;j<nclistlength(var->attributes);j++) {
		NCattribute* att = (NCattribute*)nclistget(var->attributes,j);
		ncstat = buildattribute(dapcomm,att,var->etype,varid);
        	if(ncstat != NC_NOERR) goto done;
	    }
	}
	/* Tag the variable with its DAP path */
	if(dapparamcheck(dapcomm,"show","projection"))
	    showprojection(dapcomm,var);
    }    
done:
    return THROW(ncstat);
}

static NCerror
buildglobalattrs(NCDAPCOMMON* dapcomm, CDFnode* root)
{
    int i;
    NCerror ncstat = NC_NOERR;
    const char* txt;
    char *nltxt, *p;
    NCbytes* buf = NULL;
    NClist* cdfnodes;
    NC* drno = dapcomm->controller;

    if(root->attributes != NULL) {
        for(i=0;i<nclistlength(root->attributes);i++) {
   	    NCattribute* att = (NCattribute*)nclistget(root->attributes,i);
	    ncstat = buildattribute(dapcomm,att,NC_NAT,NC_GLOBAL);
            if(ncstat != NC_NOERR) goto done;
	}
    }

    /* Add global attribute identifying the sequence dimensions */
    if(dapparamcheck(dapcomm,"show","seqdims")) {
        buf = ncbytesnew();
        cdfnodes = dapcomm->cdf.ddsroot->tree->nodes;
        for(i=0;i<nclistlength(cdfnodes);i++) {
	    CDFnode* dim = (CDFnode*)nclistget(cdfnodes,i);
	    if(dim->nctype != NC_Dimension) continue;
	    if(DIMFLAG(dim,CDFDIMSEQ)) {
	        char* cname = cdflegalname(dim->ocname);
	        if(ncbyteslength(buf) > 0) ncbytescat(buf,", ");
	        ncbytescat(buf,cname);
	        nullfree(cname);
	    }
	}
        if(ncbyteslength(buf) > 0) {
            ncstat = nc_put_att_text(drno->substrate,NC_GLOBAL,"_sequence_dimensions",
	           ncbyteslength(buf),ncbytescontents(buf));
	}
    }

    /* Define some additional system global attributes
       depending on show= clientparams*/
    /* Ignore failures*/

    if(dapparamcheck(dapcomm,"show","translate")) {
        /* Add a global attribute to show the translation */
        ncstat = nc_put_att_text(drno->substrate,NC_GLOBAL,"_translate",
	           strlen("netcdf-3"),"netcdf-3");
    }
    if(dapparamcheck(dapcomm,"show","url")) {
	if(dapcomm->oc.rawurltext != NULL)
            ncstat = nc_put_att_text(drno->substrate,NC_GLOBAL,"_url",
				       strlen(dapcomm->oc.rawurltext),dapcomm->oc.rawurltext);
    }
    if(dapparamcheck(dapcomm,"show","dds")) {
	txt = NULL;
	if(dapcomm->cdf.ddsroot != NULL)
  	    txt = oc_tree_text(dapcomm->oc.conn,dapcomm->cdf.ddsroot->ocnode);
	if(txt != NULL) {
	    /* replace newlines with spaces*/
	    nltxt = nulldup(txt);
	    for(p=nltxt;*p;p++) {if(*p == '\n' || *p == '\r' || *p == '\t') {*p = ' ';}};
            ncstat = nc_put_att_text(drno->substrate,NC_GLOBAL,"_dds",strlen(nltxt),nltxt);
	    nullfree(nltxt);
	}
    }
    if(dapparamcheck(dapcomm,"show","das")) {
	txt = NULL;
	if(dapcomm->oc.ocdasroot != NULL)
	    txt = oc_tree_text(dapcomm->oc.conn,dapcomm->oc.ocdasroot);
	if(txt != NULL) {
	    nltxt = nulldup(txt);
	    for(p=nltxt;*p;p++) {if(*p == '\n' || *p == '\r' || *p == '\t') {*p = ' ';}};
            ncstat = nc_put_att_text(drno->substrate,NC_GLOBAL,"_das",strlen(nltxt),nltxt);
	    nullfree(nltxt);
	}
    }

done:
    ncbytesfree(buf);
    return THROW(ncstat);
}

static NCerror
buildattribute(NCDAPCOMMON* dapcomm, NCattribute* att, nc_type vartype, int varid)
{
    int i;
    NCerror ncstat = NC_NOERR;
    unsigned int nvalues = nclistlength(att->values);
    NC* drno = dapcomm->controller;

    /* If the type of the attribute is string, then we need*/
    /* to convert to a single character string by concatenation.
	modified: 10/23/09 to insert newlines.
	modified: 10/28/09 to interpret escapes
    */
    if(att->etype == NC_STRING || att->etype == NC_URL) {
	char* newstring;
	size_t newlen = 0;
	for(i=0;i<nvalues;i++) {
	    char* s = (char*)nclistget(att->values,i);
	    newlen += (1+strlen(s));
	}
	newstring = (char*)malloc(newlen);
        MEMCHECK(newstring,NC_ENOMEM);
	newstring[0] = '\0';
	for(i=0;i<nvalues;i++) {
	    char* s = (char*)nclistget(att->values,i);
	    if(i > 0) strcat(newstring,"\n");
	    strcat(newstring,s);
	}
        dapexpandescapes(newstring);
	if(newstring[0]=='\0')
	    ncstat = nc_put_att_text(drno->substrate,varid,att->name,1,newstring);
	else
	    ncstat = nc_put_att_text(drno->substrate,varid,att->name,strlen(newstring),newstring);
	free(newstring);
    } else {
	nc_type atype;
	unsigned int typesize;
	void* mem;
	/* It turns out that some servers upgrade the type
           of _FillValue in order to correctly preserve the
           original value. However, since the type of the
           underlying variable is not changes, we get a type
           mismatch. So, make sure the type of the fillvalue
           is the same as that of the controlling variable.
	*/
        if(varid != NC_GLOBAL && strcmp(att->name,"_FillValue")==0)
	    atype = nctypeconvert(dapcomm,vartype);
	else
	    atype = nctypeconvert(dapcomm,att->etype);
	typesize = nctypesizeof(atype);
	mem = malloc(typesize * nvalues);
        ncstat = dapcvtattrval(atype,mem,att->values);
        ncstat = nc_put_att(drno->substrate,varid,att->name,atype,nvalues,mem);
	nullfree(mem);
    }
    return THROW(ncstat);
}

static char*
getdefinename(CDFnode* node)
{
    char* spath = NULL;
    NClist* path = NULL;

    switch (node->nctype) {
    case NC_Atomic:
	/* The define name is same as the fullname with elided nodes */
	path = nclistnew();
        collectnodepath(node,path,!WITHDATASET);
        spath = makepathstring(path,".",PATHNC|PATHELIDE);
        nclistfree(path);
	break;

    case NC_Dimension:
	/* Return just the node's ncname */
	spath = nulldup(node->ncbasename);	
	break;

    default:
	PANIC("unexpected nctype");
    }
    return spath;
}

int
NCDAP_ping(const char* url)
{
    OCerror ocstat = OC_NOERR;
    ocstat = oc_ping(url);
    return ocerrtoncerr(ocstat);
}

int
NCD2_inq_format_extended(int ncid, int* formatp, int* modep)
{
    NC* nc;
    int ncstatus = NC_check_id(ncid, (NC**)&nc);
    if(ncstatus != NC_NOERR) return THROW(ncstatus);
    if(modep) *modep = nc->mode;
    if(formatp) *formatp = NC_FORMAT_DAP2;
    return NC_NOERR;
}

/**************************************************/
/* Support functions */

/*
   Provide short and/or unified names for dimensions.
   This must mimic lib-ncdap, which is difficult.
*/
NCerror
computecdfdimnames(NCDAPCOMMON* nccomm)
{
    int i,j;
    char tmp[NC_MAX_NAME*2];
    NClist* conflicts = nclistnew();
    NClist* varnodes = nccomm->cdf.ddsroot->tree->varnodes;
    NClist* alldims;
    NClist* basedims;
    
    /* Collect all dimension nodes from dimsetall lists */

    alldims = getalldims(nccomm,0);    

    /* Assign an index to all anonymous dimensions
       vis-a-vis its containing variable
    */
    for(i=0;i<nclistlength(varnodes);i++) {
	CDFnode* var = (CDFnode*)nclistget(varnodes,i);
        for(j=0;j<nclistlength(var->array.dimsetall);j++) {
	    CDFnode* dim = (CDFnode*)nclistget(var->array.dimsetall,j);
	    if(dim->ocname != NULL) continue; /* not anonymous */
 	    computedimindexanon(dim,var);
	}
    }

    /* Unify dimensions by defining one dimension as the "base"
       dimension, and make all "equivalent" dimensions point to the
       base dimension.
	1. Equivalent means: same size and both have identical non-null names.
	2. Dims with same name but different sizes will be handled separately
    */
    for(i=0;i<nclistlength(alldims);i++) {
	CDFnode* dupdim = NULL;
	CDFnode* basedim = (CDFnode*)nclistget(alldims,i);
	if(basedim == NULL) continue;
	if(basedim->dim.basedim != NULL) continue; /* already processed*/
	for(j=i+1;j<nclistlength(alldims);j++) { /* Sigh, n**2 */
	    dupdim = (CDFnode*)nclistget(alldims,j);
	    if(basedim == dupdim) continue;
	    if(dupdim == NULL) continue;
	    if(dupdim->dim.basedim != NULL) continue; /* already processed */
	    if(!equivalentdim(basedim,dupdim))
		continue;
            dupdim->dim.basedim = basedim; /* equate */
#ifdef DEBUG1
fprintf(stderr,"assign: %s/%s -> %s/%s\n",
basedim->dim.array->ocname,basedim->ocname,
dupdim->dim.array->ocname,dupdim->ocname
);
#endif
	}
    }

    /* Next case: same name and different sizes*/
    /* => rename second dim */

    for(i=0;i<nclistlength(alldims);i++) {
	CDFnode* basedim = (CDFnode*)nclistget(alldims,i);
	if(basedim->dim.basedim != NULL) continue;
	/* Collect all conflicting dimensions */
	nclistclear(conflicts);
        for(j=i+1;j<nclistlength(alldims);j++) {
	    CDFnode* dim = (CDFnode*)nclistget(alldims,j);
	    if(dim->dim.basedim != NULL) continue;
	    if(dim->ocname == NULL && basedim->ocname == NULL) continue;
	    if(dim->ocname == NULL || basedim->ocname == NULL) continue;
	    if(strcmp(dim->ocname,basedim->ocname)!=0) continue;
	    if(dim->dim.declsize == basedim->dim.declsize) continue;
#ifdef DEBUG2
fprintf(stderr,"conflict: %s[%lu] %s[%lu]\n",
			basedim->ncfullname,(unsigned long)basedim->dim.declsize,
			dim->ncfullname,(unsigned long)dim->dim.declsize);
#endif
	    nclistpush(conflicts,(void*)dim);
	}
	/* Give  all the conflicting dimensions an index */
	for(j=0;j<nclistlength(conflicts);j++) {
	    CDFnode* dim = (CDFnode*)nclistget(conflicts,j);
	    dim->dim.index1 = j+1;
	}
    }
    nclistfree(conflicts);

    /* Replace all non-base dimensions with their base dimension */
    for(i=0;i<nclistlength(varnodes);i++) {
	CDFnode* node = (CDFnode*)nclistget(varnodes,i);
	replacedims(node->array.dimsetall);
	replacedims(node->array.dimsetplus);
	replacedims(node->array.dimset0);
    }

    /* Collect list of all basedims */
    basedims = nclistnew();
    for(i=0;i<nclistlength(alldims);i++) {
	CDFnode* dim = (CDFnode*)nclistget(alldims,i);
	if(dim->dim.basedim == NULL) {
	    if(!nclistcontains(basedims,(void*)dim)) {
		nclistpush(basedims,(void*)dim);
	    }
	}
    }

    nccomm->cdf.ddsroot->tree->dimnodes = basedims;

    /* cleanup */
    nclistfree(alldims);

    /* Assign ncbasenames and ncfullnames to base dimensions */
    for(i=0;i<nclistlength(basedims);i++) {
	CDFnode* dim = (CDFnode*)nclistget(basedims,i);
	CDFnode* var = dim->dim.array;
	if(dim->dim.basedim != NULL) PANIC1("nonbase basedim: %s\n",dim->ocname);
	/* stringdim names are already assigned */
	if(dim->ocname == NULL) { /* anonymous: use the index to compute the name */
            snprintf(tmp,sizeof(tmp),"%s_%d",
                            var->ncfullname,dim->dim.index1-1);
            nullfree(dim->ncbasename);
            dim->ncbasename = cdflegalname(tmp);
            nullfree(dim->ncfullname);
            dim->ncfullname = nulldup(dim->ncbasename);
    	} else { /* !anonymous; use index1 if defined */
   	    char* legalname = cdflegalname(dim->ocname);
	    nullfree(dim->ncbasename);
	    if(dim->dim.index1 > 0) {/* need to fix conflicting names (see above) */
	        char sindex[64];
		snprintf(sindex,sizeof(sindex),"_%d",dim->dim.index1);
		dim->ncbasename = (char*)malloc(strlen(sindex)+strlen(legalname)+1);
		if(dim->ncbasename == NULL) {nullfree(legalname); return NC_ENOMEM;}
		strcpy(dim->ncbasename,legalname);
		strcat(dim->ncbasename,sindex);
		nullfree(legalname);
	    } else {/* standard case */
	        dim->ncbasename = legalname;
	    }
    	    nullfree(dim->ncfullname);
	    dim->ncfullname = nulldup(dim->ncbasename);
	}
     }

    /* Verify unique and defined names for dimensions*/
    for(i=0;i<nclistlength(basedims);i++) {
	CDFnode* dim1 = (CDFnode*)nclistget(basedims,i);
	if(dim1->dim.basedim != NULL) PANIC1("nonbase basedim: %s\n",dim1->ncbasename);
	if(dim1->ncbasename == NULL || dim1->ncfullname == NULL)
	    PANIC1("missing dim names: %s",dim1->ocname);
	/* search backward so we can delete duplicates */
	for(j=nclistlength(basedims)-1;j>i;j--) {
	    CDFnode* dim2 = (CDFnode*)nclistget(basedims,j);
	    if(strcmp(dim1->ncfullname,dim2->ncfullname)==0) {
		/* complain and suppress one of them */
		fprintf(stderr,"duplicate dim names: %s[%lu] %s[%lu]\n",
			dim1->ncfullname,(unsigned long)dim1->dim.declsize,
			dim2->ncfullname,(unsigned long)dim2->dim.declsize);
		nclistremove(basedims,j);
	    }
	}
    }

#ifdef DEBUG
for(i=0;i<nclistlength(basedims);i++) {
CDFnode* dim = (CDFnode*)nclistget(basedims,i);
fprintf(stderr,"basedim: %s=%ld\n",dim->ncfullname,(long)dim->dim.declsize);
 }
#endif

    return NC_NOERR;
}

int
constrainable(NCURI* durl)
{
   char** protocol = constrainableprotocols;
   for(;*protocol;protocol++) {
	if(strcmp(durl->protocol,*protocol)==0)
	    return 1;
   }
   return 0;
}

/* Note: this routine only applies some common
   client parameters, other routines may apply
   specific ones.
*/

static NCerror
applyclientparams(NCDAPCOMMON* nccomm)
{
    int i,len;
    int dfaltstrlen = DEFAULTSTRINGLENGTH;
    int dfaltseqlim = DEFAULTSEQLIMIT;
    const char* value;
    char tmpname[NC_MAX_NAME+32];
    char* pathstr;
    OClink conn = nccomm->oc.conn;
    unsigned long limit;

    ASSERT(nccomm->oc.url != NULL);

    nccomm->cdf.cache->cachelimit = DFALTCACHELIMIT;
    value = oc_clientparam_get(conn,"cachelimit");
    limit = getlimitnumber(value);
    if(limit > 0) nccomm->cdf.cache->cachelimit = limit;

    nccomm->cdf.fetchlimit = DFALTFETCHLIMIT;
    value = oc_clientparam_get(conn,"fetchlimit");
    limit = getlimitnumber(value);
    if(limit > 0) nccomm->cdf.fetchlimit = limit;

    nccomm->cdf.smallsizelimit = DFALTSMALLLIMIT;
    value = oc_clientparam_get(conn,"smallsizelimit");
    limit = getlimitnumber(value);
    if(limit > 0) nccomm->cdf.smallsizelimit = limit;

    nccomm->cdf.cache->cachecount = DFALTCACHECOUNT;
#ifdef HAVE_GETRLIMIT
    { struct rlimit rl;
      if(getrlimit(RLIMIT_NOFILE, &rl) >= 0) {
	nccomm->cdf.cache->cachecount = (size_t)(rl.rlim_cur / 2);
      }
    }
#endif
    value = oc_clientparam_get(conn,"cachecount");
    limit = getlimitnumber(value);
    if(limit > 0) nccomm->cdf.cache->cachecount = limit;
    /* Ignore limit if not caching */
    if(!FLAGSET(nccomm->controls,NCF_CACHE))
        nccomm->cdf.cache->cachecount = 0;

    if(oc_clientparam_get(conn,"nolimit") != NULL)
	dfaltseqlim = 0;
    value = oc_clientparam_get(conn,"limit");
    if(value != NULL && strlen(value) != 0) {
        if(sscanf(value,"%d",&len) && len > 0) dfaltseqlim = len;
    }
    nccomm->cdf.defaultsequencelimit = dfaltseqlim;

    /* allow embedded _ */
    value = oc_clientparam_get(conn,"stringlength");
    if(value != NULL && strlen(value) != 0) {
        if(sscanf(value,"%d",&len) && len > 0) dfaltstrlen = len;
    }
    nccomm->cdf.defaultstringlength = dfaltstrlen;

    /* String dimension limits apply to variables */
    for(i=0;i<nclistlength(nccomm->cdf.ddsroot->tree->varnodes);i++) {
	CDFnode* var = (CDFnode*)nclistget(nccomm->cdf.ddsroot->tree->varnodes,i);
	/* Define the client param stringlength for this variable*/
	var->maxstringlength = 0; /* => use global dfalt */
	strcpy(tmpname,"stringlength_");
	pathstr = makeocpathstring(conn,var->ocnode,".");
	strncat(tmpname,pathstr,NC_MAX_NAME);
	nullfree(pathstr);
	value = oc_clientparam_get(conn,tmpname);	
        if(value != NULL && strlen(value) != 0) {
            if(sscanf(value,"%d",&len) && len > 0) var->maxstringlength = len;
	}
    }
    /* Sequence limits apply to sequences */
    for(i=0;i<nclistlength(nccomm->cdf.ddsroot->tree->nodes);i++) {
	CDFnode* var = (CDFnode*)nclistget(nccomm->cdf.ddsroot->tree->nodes,i);
	if(var->nctype != NC_Sequence) continue;
	var->sequencelimit = dfaltseqlim;
	strcpy(tmpname,"nolimit_");
	pathstr = makeocpathstring(conn,var->ocnode,".");
	strncat(tmpname,pathstr,NC_MAX_NAME);
	if(oc_clientparam_get(conn,tmpname) != NULL)
	    var->sequencelimit = 0;
	strcpy(tmpname,"limit_");
	strncat(tmpname,pathstr,NC_MAX_NAME);
	value = oc_clientparam_get(conn,tmpname);
        if(value != NULL && strlen(value) != 0) {
            if(sscanf(value,"%d",&len) && len > 0)
		var->sequencelimit = len;
	}
	nullfree(pathstr);
    }

    /* test for the appropriate fetch flags */
    value = oc_clientparam_get(conn,"fetch");
    if(value != NULL && strlen(value) > 0) {
	if(value[0] == 'd' || value[0] == 'D') {
            SETFLAG(nccomm->controls,NCF_ONDISK);
	}
    }

    /* test for the force-whole-var flag */
    value = oc_clientparam_get(conn,"wholevar");
    if(value != NULL) {
        SETFLAG(nccomm->controls,NCF_WHOLEVAR);
    }

    return NC_NOERR;
}

/**
 *  Given an anonymous dimension, compute the
 *  effective 0-based index wrt to the specified var.
 *  The result should mimic the libnc-dap indices.
 */

static void
computedimindexanon(CDFnode* dim, CDFnode* var)
{
    int i;
    NClist* dimset = var->array.dimsetall;
    for(i=0;i<nclistlength(dimset);i++) {
	CDFnode* candidate = (CDFnode*)nclistget(dimset,i);
        if(dim == candidate) {
	   dim->dim.index1=i+1;
	   return;
	}
    }
}

/* Replace dims in a list with their corresponding basedim */
static void
replacedims(NClist* dims)
{
    int i;
    for(i=0;i<nclistlength(dims);i++) {
        CDFnode* dim = (CDFnode*)nclistget(dims,i);
	CDFnode* basedim = dim->dim.basedim;
	if(basedim == NULL) continue;
	nclistset(dims,i,(void*)basedim);
    }
}

/**
 Two dimensions are equivalent if
 1. they have the same size
 2. neither are anonymous
 3. they ave the same names. 
 */
static int
equivalentdim(CDFnode* basedim, CDFnode* dupdim)
{
    if(dupdim->dim.declsize != basedim->dim.declsize) return 0;
    if(basedim->ocname == NULL && dupdim->ocname == NULL) return 0;
    if(basedim->ocname == NULL || dupdim->ocname == NULL) return 0;
    if(strcmp(dupdim->ocname,basedim->ocname) != 0) return 0;
    return 1;
}

static void
getalldimsa(NClist* dimset, NClist* alldims)
{
    int i;
    for(i=0;i<nclistlength(dimset);i++) {
	CDFnode* dim = (CDFnode*)nclistget(dimset,i);
	if(!nclistcontains(alldims,(void*)dim)) {
#ifdef DEBUG3
fprintf(stderr,"getalldims: %s[%lu]\n",
			dim->ncfullname,(unsigned long)dim->dim.declsize);
#endif
	    nclistpush(alldims,(void*)dim);
	}
    }
}

/* Accumulate a set of all the known dimensions
   vis-a-vis defined variables
*/
NClist*
getalldims(NCDAPCOMMON* nccomm, int visibleonly)
{
    int i;
    NClist* alldims = nclistnew();
    NClist* varnodes = nccomm->cdf.ddsroot->tree->varnodes;

    /* get bag of all dimensions */
    for(i=0;i<nclistlength(varnodes);i++) {
	CDFnode* node = (CDFnode*)nclistget(varnodes,i);
	if(!visibleonly || !node->invisible) {
	    getalldimsa(node->array.dimsetall,alldims);
	}
    }
    return alldims;
}


static NCerror
addstringdims(NCDAPCOMMON* dapcomm)
{
    /* for all variables of string type, we will need another dimension
       to represent the string; Accumulate the needed sizes and create
       the dimensions with a specific name: either as specified
       in DODS{...} attribute set or defaulting to the variable name.
       All such dimensions are global.
    */
    int i;
    NClist* varnodes = dapcomm->cdf.ddsroot->tree->varnodes;
    CDFnode* globalsdim = NULL;
    char dimname[4096];
    size_t dimsize;

    /* Start by creating the global string dimension */
    snprintf(dimname,sizeof(dimname),"maxStrlen%lu",
	    (unsigned long)dapcomm->cdf.defaultstringlength);
    globalsdim = makecdfnode(dapcomm, dimname, OC_Dimension, NULL,
                                 dapcomm->cdf.ddsroot);
    nclistpush(dapcomm->cdf.ddsroot->tree->nodes,(void*)globalsdim);
    DIMFLAGSET(globalsdim,CDFDIMSTRING);
    globalsdim->dim.declsize = dapcomm->cdf.defaultstringlength;
    globalsdim->dim.declsize0 = globalsdim->dim.declsize;
    globalsdim->dim.array = dapcomm->cdf.ddsroot;
    globalsdim->ncbasename = cdflegalname(dimname);
    globalsdim->ncfullname = nulldup(globalsdim->ncbasename);
    dapcomm->cdf.globalstringdim = globalsdim;

    for(i=0;i<nclistlength(varnodes);i++) {
	CDFnode* var = (CDFnode*)nclistget(varnodes,i);
	CDFnode* sdim = NULL;

	/* Does this node need a string dim? */
	if(var->etype != NC_STRING && var->etype != NC_URL) continue;

	dimsize = 0;
	if(var->dodsspecial.maxstrlen > 0)
	    dimsize = var->dodsspecial.maxstrlen;
	else
	    dimsize = var->maxstringlength;

	/* check is a variable-specific string length was specified */
	if(dimsize == 0) 
	    sdim = dapcomm->cdf.globalstringdim; /* use default */
	else {
	    /* create a psuedo dimension for the charification of the string*/
	    if(var->dodsspecial.dimname != NULL) {
	        strncpy(dimname,var->dodsspecial.dimname,sizeof(dimname));
	        dimname[sizeof(dimname)-1] = '\0';
	    } else
	        snprintf(dimname,sizeof(dimname),"maxStrlen%lu",
			 (unsigned long)dimsize);
	    sdim = makecdfnode(dapcomm, dimname, OC_Dimension, NULL,
                                 dapcomm->cdf.ddsroot);
	    if(sdim == NULL) return THROW(NC_ENOMEM);
	    nclistpush(dapcomm->cdf.ddsroot->tree->nodes,(void*)sdim);
	    DIMFLAGSET(sdim,CDFDIMSTRING);
	    sdim->dim.declsize = dimsize;
	    sdim->dim.declsize0 = dimsize;
	    sdim->dim.array = var;
	    sdim->ncbasename = cdflegalname(sdim->ocname);
	    sdim->ncfullname = nulldup(sdim->ncbasename);
	}
	/* tag the variable with its string dimension*/
	var->array.stringdim = sdim;
    }
    return NC_NOERR;
}

static NCerror
defrecorddim(NCDAPCOMMON* dapcomm)
{
    unsigned int i;
    NCerror ncstat = NC_NOERR;
    NClist* basedims;

    if(dapcomm->cdf.recorddimname == NULL) return NC_NOERR; /* ignore */
    /* Locate the base dimension matching the record dim */
    basedims = dapcomm->cdf.ddsroot->tree->dimnodes;
    for(i=0;i<nclistlength(basedims);i++) {
        CDFnode* dim = (CDFnode*)nclistget(basedims,i);
	if(strcmp(dim->ocname,dapcomm->cdf.recorddimname) != 0) continue;
        DIMFLAGSET(dim,CDFDIMRECORD);
	dapcomm->cdf.recorddim = dim;
	break;
    }

    return ncstat;
}

static NCerror
defseqdims(NCDAPCOMMON* dapcomm)
{
    unsigned int i = 0;
    NCerror ncstat = NC_NOERR;
    int seqdims = 1; /* default is to compute seq dims counts */

    /* Does the user want to compute actual sequence sizes? */
    if(dapparamvalue(dapcomm,"noseqdims")) seqdims = 0;

    /*
	Compute and define pseudo dimensions for sequences
	meeting the following qualifications:
	1. all parents (transitively) of the sequence must
           be either a dataset or a scalar structure.
	2. it must be possible to find a usable sequence constraint.
	All other sequences will be ignored.
    */

    for(i=0;i<nclistlength(dapcomm->cdf.ddsroot->tree->seqnodes);i++) {
        CDFnode* seq = (CDFnode*)nclistget(dapcomm->cdf.ddsroot->tree->seqnodes,i);
	size_t seqsize = 0;
	CDFnode* sqdim = NULL;
	CDFnode* container;
	/* Does this sequence match the requirements for use ? */
	seq->usesequence = 1; /* assume */
	for(container=seq->container;container != NULL;container=container->container) {
	    if(container->nctype == NC_Dataset) break;
	    if(container->nctype != NC_Structure
	       || nclistlength(container->array.dimset0) > 0)
		{seq->usesequence = 0; break;}/* no good */
	}	
	/* Does the user want us to compute the actual sequence dim size? */
	if(seq->usesequence && seqdims) {
	    ncstat = getseqdimsize(dapcomm,seq,&seqsize);
	    if(ncstat != NC_NOERR) {
                /* Cannot read sequence; mark as unusable */
	        seq->usesequence = 0;
	    }
	} else { /* !seqdims default to size = 1 */
	    seqsize = 1; 
	}
	if(seq->usesequence) {
	    /* Note: we are making the dimension in the dds root tree */
            ncstat = makeseqdim(dapcomm,seq,seqsize,&sqdim);
            if(ncstat) goto fail;
            seq->array.seqdim = sqdim;
	} else
            seq->array.seqdim = NULL;
    }

fail:
    return ncstat;
}

static NCerror
showprojection(NCDAPCOMMON* dapcomm, CDFnode* var)
{
    int i,rank;
    NCerror ncstat = NC_NOERR;
    NCbytes* projection = ncbytesnew();
    NClist* path = nclistnew();
    NC* drno = dapcomm->controller;

    /* Collect the set of DDS node name forming the xpath */
    collectnodepath(var,path,WITHOUTDATASET);
    for(i=0;i<nclistlength(path);i++) {
        CDFnode* node = (CDFnode*)nclistget(path,i);
	if(i > 0) ncbytescat(projection,".");
	ncbytescat(projection,node->ocname);
    }
    /* Now, add the dimension info */
    rank = nclistlength(var->array.dimset0);
    for(i=0;i<rank;i++) {
	CDFnode* dim = (CDFnode*)nclistget(var->array.dimset0,i);
	char tmp[32];
	ncbytescat(projection,"[");
	snprintf(tmp,sizeof(tmp),"%lu",(unsigned long)dim->dim.declsize);
	ncbytescat(projection,tmp);
	ncbytescat(projection,"]");
    }    
    /* Define the attribute */
    ncstat = nc_put_att_text(getncid(drno),var->ncid,
                               "_projection",
		               ncbyteslength(projection),
			       ncbytescontents(projection));
    return ncstat;
}

static NCerror
getseqdimsize(NCDAPCOMMON* dapcomm, CDFnode* seq, size_t* sizep)
{
    NCerror ncstat = NC_NOERR;
    OCerror ocstat = OC_NOERR;
    OClink conn = dapcomm->oc.conn;
    OCdatanode rootcontent = NULL;
    OCddsnode ocroot;
    CDFnode* dxdroot;
    CDFnode* xseq;
    NCbytes* seqcountconstraints = ncbytesnew();
    size_t seqsize = 0;

    /* Read the minimal amount of data in order to get the count */
    /* If the url is unconstrainable, then get the whole thing */
    computeseqcountconstraints(dapcomm,seq,seqcountconstraints);
#ifdef DEBUG
fprintf(stderr,"seqcountconstraints: %s\n",ncbytescontents(seqcountconstraints));
#endif

    /* Fetch the minimal data */
    if(FLAGSET(dapcomm->controls,NCF_UNCONSTRAINABLE))
        ncstat = dap_fetch(dapcomm,conn,NULL,OCDATADDS,&ocroot);
    else
        ncstat = dap_fetch(dapcomm,conn,ncbytescontents(seqcountconstraints),OCDATADDS,&ocroot);
    if(ncstat) goto fail;

    ncstat = buildcdftree(dapcomm,ocroot,OCDATA,&dxdroot);
    if(ncstat) goto fail;	
    /* attach DATADDS to DDS */
    ncstat = attach(dxdroot,seq);
    if(ncstat) goto fail;	

    /* WARNING: we are now switching to datadds tree */
    xseq = seq->attachment;
    ncstat = countsequence(dapcomm,xseq,&seqsize);
    if(ncstat) goto fail;

#ifdef DEBUG
fprintf(stderr,"sequencesize: %s = %lu\n",seq->ocname,(unsigned long)seqsize);
#endif

    /* throw away the fetch'd trees */
    unattach(dapcomm->cdf.ddsroot);
    freecdfroot(dxdroot);
    if(ncstat != NC_NOERR) {
        /* Cannot get DATADDDS*/
	char* code;
	char* msg;
	long httperr;
	oc_svcerrordata(dapcomm->oc.conn,&code,&msg,&httperr);
	if(code != NULL) {
	    nclog(NCLOGERR,"oc_fetch_datadds failed: %s %s %l",
			code,msg,httperr);
	}
	ocstat = OC_NOERR;
    }		
    if(sizep) *sizep = seqsize;

fail:
    ncbytesfree(seqcountconstraints);
    oc_data_free(conn,rootcontent);
    if(ocstat) ncstat = ocerrtoncerr(ocstat);
    return ncstat;
}

static NCerror
makeseqdim(NCDAPCOMMON* dapcomm, CDFnode* seq, size_t count, CDFnode** sqdimp)
{
    CDFnode* sqdim;
    CDFnode* root = seq->root;
    CDFtree* tree = root->tree;

    /* build the dimension with given size; keep the dimension anonymous */
    sqdim = makecdfnode(dapcomm,seq->ocname,OC_Dimension,NULL,root);
    if(sqdim == NULL) return THROW(NC_ENOMEM);
    nclistpush(tree->nodes,(void*)sqdim);
    /* Assign a name to the sequence node */
    sqdim->ncbasename = cdflegalname(seq->ocname);
    sqdim->ncfullname = nulldup(sqdim->ncbasename);
    DIMFLAGSET(sqdim,CDFDIMSEQ);
    sqdim->dim.declsize = count;
    sqdim->dim.declsize0 = count;
    sqdim->dim.array = seq;
    if(sqdimp) *sqdimp = sqdim;
    return NC_NOERR;
}

static NCerror
countsequence(NCDAPCOMMON* dapcomm, CDFnode* xseq, size_t* sizep)
{
    unsigned int i;
    NClist* path = nclistnew();
    int index;
    OCerror ocstat = OC_NOERR;
    NCerror ncstat = NC_NOERR;
    OClink conn = dapcomm->oc.conn;
    size_t recordcount;
    CDFnode* xroot;
    OCdatanode data = NULL;

    ASSERT((xseq->nctype == NC_Sequence));

    /* collect the path to the sequence node */
    collectnodepath(xseq,path,WITHDATASET);

    /* Get tree root */
    ASSERT(xseq->root == (CDFnode*)nclistget(path,0));
    xroot = xseq->root;
    ocstat = oc_data_getroot(conn,xroot->tree->ocroot,&data);
    if(ocstat) goto done;

    /* Basically we use the path to walk the data instances to reach
       the sequence instance
    */
    for(i=0;i<nclistlength(path);i++) {
        CDFnode* current = (CDFnode*)nclistget(path,i);
	OCdatanode nextdata = NULL;
	CDFnode* next = NULL;
	
	/* invariant: current = ith node in path; data = corresponding
           datanode
        */

	/* get next node in next and next instance in nextdata */
	if(current->nctype == NC_Structure
	   || current->nctype == NC_Dataset) {
	    if(nclistlength(current->array.dimset0) > 0) {
		/* Cannot handle this case */
		ncstat = THROW(NC_EDDS);
		goto done;
	    }
	    /* get next node in path; structure/dataset => exists */
	    next = (CDFnode*)nclistget(path,i+1);
	    index = fieldindex(current,next);
            /* Move to appropriate field */
	    ocstat = oc_data_ithfield(conn,data,index,&nextdata);
	    if(ocstat) goto done;
	    oc_data_free(conn,data);
	    data = nextdata; /* set up for next loop iteration */
	} else if(current->nctype ==  NC_Sequence) {
	    /* Check for nested Sequences */
	    if(current != xseq) {
		/* Cannot handle this case */
		ncstat = THROW(NC_EDDS);
		goto done;
	    }
	    /* Get the record count */
	    ocstat = oc_data_recordcount(conn,data,&recordcount);
    	    if(sizep) *sizep = recordcount;
	    oc_data_free(conn,data); /* reclaim */
	    break; /* leave the loop */
	} else {
	    PANIC("unexpected mode");
	    return NC_EINVAL;
        }
    }

done:
    nclistfree(path);
    if(ocstat) ncstat = ocerrtoncerr(ocstat);
    return THROW(ncstat);
}

static NCerror
freeNCDAPCOMMON(NCDAPCOMMON* dapcomm)
{
    freenccache(dapcomm,dapcomm->cdf.cache);
    nclistfree(dapcomm->cdf.projectedvars);
    nullfree(dapcomm->cdf.recorddimname);

    /* free the trees */
    freecdfroot(dapcomm->cdf.ddsroot);
    dapcomm->cdf.ddsroot = NULL;
    freecdfroot(dapcomm->cdf.fullddsroot);
    dapcomm->cdf.fullddsroot = NULL;
    if(dapcomm->oc.ocdasroot != NULL)
        oc_root_free(dapcomm->oc.conn,dapcomm->oc.ocdasroot);
    dapcomm->oc.ocdasroot = NULL;
    oc_close(dapcomm->oc.conn); /* also reclaims remaining OC trees */
    ncurifree(dapcomm->oc.url);
    nullfree(dapcomm->oc.urltext);
    nullfree(dapcomm->oc.rawurltext);

    dcefree((DCEnode*)dapcomm->oc.dapconstraint);
    dapcomm->oc.dapconstraint = NULL;

    free(dapcomm);

    return NC_NOERR;
}

static int
fieldindex(CDFnode* parent, CDFnode* child)
{
    unsigned int i;
    for(i=0;i<nclistlength(parent->subnodes);i++) {
	CDFnode* node = (CDFnode*)nclistget(parent->subnodes,i);
	if(node == child) return i;
    }
    return -1;
}

/*
This is more complex than one might think. We want to find
a path to a variable inside the given node so that we can
ask for a single instance of that variable to minimize the
amount of data we retrieve. However, we want to avoid passing
through any nested sequence. This is possible because of the way
that sequencecheck() works.
TODO: some servers will not accept an unconstrained fetch, so
make sure we always have a constraint.
*/

static NCerror
computeseqcountconstraints(NCDAPCOMMON* dapcomm, CDFnode* seq, NCbytes* seqcountconstraints)
{
    int i,j;
    NClist* path = NULL;
    CDFnode* var = NULL;

    ASSERT(seq->nctype == NC_Sequence);
    computeseqcountconstraintsr(dapcomm,seq,&var);

    ASSERT((var != NULL));

    /* Compute var path */
    path = nclistnew();
    collectnodepath(var,path,WITHOUTDATASET);

    /* construct the projection path using minimal index values */
    for(i=0;i<nclistlength(path);i++) {
	CDFnode* node = (CDFnode*)nclistget(path,i);
	if(i > 0) ncbytescat(seqcountconstraints,".");
	ncbytescat(seqcountconstraints,node->ocname);
	if(node == seq) {
	    /* Use the limit */
	    if(node->sequencelimit > 0) {
		char tmp[64];
		snprintf(tmp,sizeof(tmp),"[0:%lu]",
		         (unsigned long)(node->sequencelimit - 1));
		ncbytescat(seqcountconstraints,tmp);
	    }
	} else if(nclistlength(node->array.dimset0) > 0) {
	    int ndims = nclistlength(node->array.dimset0);
	    for(j=0;j<ndims;j++) {
		CDFnode* dim = (CDFnode*)nclistget(node->array.dimset0,j);
		if(DIMFLAG(dim,CDFDIMSTRING)) {
		    ASSERT((j == (ndims - 1)));
		    break;
		}
		ncbytescat(seqcountconstraints,"[0]");
	    }
	}
    }
    /* Finally, add in any selection from the original URL */
    if(dapcomm->oc.url->selection != NULL)
        ncbytescat(seqcountconstraints,dapcomm->oc.url->selection);
    nclistfree(path);
    return NC_NOERR;    
}


/* Given an existing candidate, see if we prefer newchoice */
static CDFnode*
prefer(CDFnode* candidate, CDFnode* newchoice)
{
    nc_type newtyp;
    nc_type cantyp;
    int newisstring;
    int canisstring;
    int newisscalar;
    int canisscalar;

    /* always choose !null over null */
    if(newchoice == NULL)
	return candidate;
    if(candidate == NULL)
	return newchoice;

    newtyp = newchoice->etype;
    cantyp = candidate->etype;
    newisstring = (newtyp == NC_STRING || newtyp == NC_URL);
    canisstring = (cantyp == NC_STRING || cantyp == NC_URL);
    newisscalar = (nclistlength(newchoice->array.dimset0) == 0);
    canisscalar = (nclistlength(candidate->array.dimset0) == 0);

    ASSERT(candidate->nctype == NC_Atomic && newchoice->nctype == NC_Atomic);
    
    /* choose non-string over string */
    if(canisstring && !newisstring)
	return newchoice;
    if(!canisstring && newisstring)
	return candidate;

    /* choose scalar over array */
    if(canisscalar && !newisscalar)
	return candidate;
    if(!canisscalar && newisscalar)
	return candidate;

    /* otherwise choose existing candidate */
    return candidate;
}

/* computeseqcountconstraints recursive helper function */
static void
computeseqcountconstraintsr(NCDAPCOMMON* dapcomm, CDFnode* node, CDFnode** candidatep)
{
    CDFnode* candidate;
    CDFnode* compound;
    unsigned int i;

    candidate = NULL;
    compound = NULL;
    if(node == NULL)
      return;
    for(i=0;i<nclistlength(node->subnodes);i++) {
        CDFnode* subnode = (CDFnode*)nclistget(node->subnodes,i);
        if(subnode->nctype == NC_Structure || subnode->nctype == NC_Grid)
	    compound = subnode; /* save for later recursion */
	else if(subnode->nctype == NC_Atomic)
	    candidate = prefer(candidate,subnode);
    }
    if(candidate == NULL && compound == NULL) {
	PANIC("cannot find candidate for seqcountconstraints for a sequence");
    } else if(candidate != NULL && candidatep != NULL) {
	*candidatep = candidate;
    } else { /* compound != NULL by construction */
	/* recurse on a nested grids or strucures */
        computeseqcountconstraintsr(dapcomm,compound,candidatep);
    }
}


static unsigned long
cdftotalsize(NClist* dimensions)
{
    unsigned int i;
    unsigned long total = 1;
    if(dimensions != NULL) {
	for(i=0;i<nclistlength(dimensions);i++) {
	    CDFnode* dim = (CDFnode*)nclistget(dimensions,i);
	    total *= dim->dim.declsize;
	}
    }
    return total;
}

/* Estimate variables sizes and then resort the variable list
   by that size
*/
static void
estimatevarsizes(NCDAPCOMMON* dapcomm)
{
    int ivar;
    unsigned int rank;
    size_t totalsize = 0;

    for(ivar=0;ivar<nclistlength(dapcomm->cdf.ddsroot->tree->varnodes);ivar++) {
        CDFnode* var = (CDFnode*)nclistget(dapcomm->cdf.ddsroot->tree->varnodes,ivar);
	NClist* ncdims = var->array.dimset0;
	rank = nclistlength(ncdims);
	if(rank == 0) { /* use instance size of the type */
	    var->estimatedsize = nctypesizeof(var->etype);
#ifdef DEBUG1
fprintf(stderr,"scalar %s.estimatedsize = %lu\n",
	makecdfpathstring(var,"."),var->estimatedsize);
#endif
	} else {
	    unsigned long size = cdftotalsize(ncdims);
	    size *= nctypesizeof(var->etype);
#ifdef DEBUG1
fprintf(stderr,"array %s(%u).estimatedsize = %lu\n",
	makecdfpathstring(var,"."),rank,size);
#endif
	    var->estimatedsize = size;
	}
	totalsize += var->estimatedsize;
    }
#ifdef DEBUG1
fprintf(stderr,"total estimatedsize = %lu\n",totalsize);
#endif
    dapcomm->cdf.totalestimatedsize = totalsize;
}

static NCerror
fetchtemplatemetadata(NCDAPCOMMON* dapcomm)
{
    NCerror ncstat = NC_NOERR;
    OCerror ocstat = OC_NOERR;
    OCddsnode ocroot = NULL;
    CDFnode* ddsroot = NULL;
    char* ce = NULL;

    /* Temporary hack: we need to get the selection string
       from the url
    */
    /* Get (almost) unconstrained DDS; In order to handle functions
       correctly, those selections must always be included
    */
    if(FLAGSET(dapcomm->controls,NCF_UNCONSTRAINABLE))
	ce = NULL;
    else
        ce = nulldup(dapcomm->oc.url->selection);

    /* Get selection constrained DDS */
    ncstat = dap_fetch(dapcomm,dapcomm->oc.conn,ce,OCDDS,&ocroot);
    if(ncstat != NC_NOERR) {
	/* Special Hack. If the protocol is file, then see if
           we can get the dds from the .dods file
        */
	if(strcmp(dapcomm->oc.url->protocol,"file") != 0) {
	    THROWCHK(ocstat); goto done;
	}
	/* Fetch the data dds */
        ncstat = dap_fetch(dapcomm,dapcomm->oc.conn,ce,OCDATADDS,&ocroot);
        if(ncstat != NC_NOERR) {
	    THROWCHK(ncstat); goto done;
	}
	/* Note what we did */
	nclog(NCLOGWARN,"Cannot locate .dds file, using .dods file");
    }

    /* Get selection constrained DAS */
    ncstat = dap_fetch(dapcomm,dapcomm->oc.conn,ce,OCDAS,&dapcomm->oc.ocdasroot);
    if(ncstat != NC_NOERR) {
	/* Ignore but complain */
	nclog(NCLOGWARN,"Could not read DAS; ignored");
        dapcomm->oc.ocdasroot = NULL;	
	ncstat = NC_NOERR;
    }

    /* Construct the netcdf cdf tree corresponding to the dds tree*/
    ncstat = buildcdftree(dapcomm,ocroot,OCDDS,&ddsroot);
    if(ncstat != NC_NOERR) {THROWCHK(ncstat); goto done;}
    dapcomm->cdf.fullddsroot = ddsroot;
    ddsroot = NULL; /* avoid double reclaim */

    /* Combine DDS and DAS */
    if(dapcomm->oc.ocdasroot != NULL) {
	ncstat = dapmerge(dapcomm,dapcomm->cdf.fullddsroot,
                           dapcomm->oc.ocdasroot);
        if(ncstat != NC_NOERR) {THROWCHK(ncstat); goto done;}
    }

#ifdef DEBUG
fprintf(stderr,"full template:\n%s",dumptree(dapcomm->cdf.fullddsroot));
#endif

done:
    nullfree(ce);
    if(ocstat != OC_NOERR) ncstat = ocerrtoncerr(ocstat);
    return ncstat;
}

static NCerror
fetchconstrainedmetadata(NCDAPCOMMON* dapcomm)
{
    NCerror ncstat = NC_NOERR;
    OCerror ocstat = OC_NOERR;
    OCddsnode ocroot;
    CDFnode* ddsroot; /* constrained */
    char* ce = NULL;

    if(FLAGSET(dapcomm->controls,NCF_UNCONSTRAINABLE))
	ce = NULL;
    else
        ce = buildconstraintstring(dapcomm->oc.dapconstraint);
    {
        ncstat = dap_fetch(dapcomm,dapcomm->oc.conn,ce,OCDDS,&ocroot);
        if(ncstat != NC_NOERR) {THROWCHK(ncstat); goto fail;}

        /* Construct our parallel dds tree; including attributes*/
        ncstat = buildcdftree(dapcomm,ocroot,OCDDS,&ddsroot);
        if(ncstat) goto fail;
	ocroot = NULL; /* avoid duplicate reclaim */

	dapcomm->cdf.ddsroot = ddsroot;
	ddsroot = NULL; /* to avoid double reclamation */

        if(!FLAGSET(dapcomm->controls,NCF_UNCONSTRAINABLE)) {
            /* fix DAP server problem by adding back any inserting needed structure nodes */
            ncstat = restruct(dapcomm, dapcomm->cdf.ddsroot,dapcomm->cdf.fullddsroot,dapcomm->oc.dapconstraint->projections);    
            if(ncstat) goto fail;
	}

#ifdef DEBUG
fprintf(stderr,"constrained:\n%s",dumptree(dapcomm->cdf.ddsroot));
#endif

        /* Combine DDS and DAS */
	if(dapcomm->oc.ocdasroot != NULL) {
            ncstat = dapmerge(dapcomm,dapcomm->cdf.ddsroot,
                               dapcomm->oc.ocdasroot);
            if(ncstat != NC_NOERR) {THROWCHK(ncstat); goto fail;}
	}

        /* map the constrained DDS to the unconstrained DDS */
        ncstat = mapnodes(dapcomm->cdf.ddsroot,dapcomm->cdf.fullddsroot);
        if(ncstat) goto fail;

    }

fail:
    nullfree(ce);
    if(ocstat != OC_NOERR) ncstat = ocerrtoncerr(ocstat);
    return ncstat;
}

/* Suppress variables not in usable sequences*/
static NCerror
suppressunusablevars(NCDAPCOMMON* dapcomm)
{
    int i,j;
    int found = 1;
    NClist* path = nclistnew();

    while(found) {
	found = 0;
	/* Walk backwards to aid removal semantics */
	for(i=nclistlength(dapcomm->cdf.ddsroot->tree->varnodes)-1;i>=0;i--) {
	    CDFnode* var = (CDFnode*)nclistget(dapcomm->cdf.ddsroot->tree->varnodes,i);
	    /* See if this var is under an unusable sequence */
	    nclistclear(path);
	    collectnodepath(var,path,WITHOUTDATASET);
	    for(j=0;j<nclistlength(path);j++) {
		CDFnode* node = (CDFnode*)nclistget(path,j);
		if(node->nctype == NC_Sequence
		   && !node->usesequence) {
#ifdef DEBUG
fprintf(stderr,"suppressing var in unusable sequence: %s.%s\n",node->ncfullname,var->ncbasename);
#endif
		    found = 1;
		    break;
		}
	    }
	    if(found) break;
	}
        if(found) nclistremove(dapcomm->cdf.ddsroot->tree->varnodes,i);
    }
    nclistfree(path);
    return NC_NOERR;
}


/*
For variables which have a zero size dimension,
make them invisible.
*/
static NCerror
fixzerodims(NCDAPCOMMON* dapcomm)
{
    int i,j;
    for(i=0;i<nclistlength(dapcomm->cdf.ddsroot->tree->varnodes);i++) {
	CDFnode* var = (CDFnode*)nclistget(dapcomm->cdf.ddsroot->tree->varnodes,i);
        NClist* ncdims = var->array.dimsetplus;
	if(nclistlength(ncdims) == 0) continue;
        for(j=0;j<nclistlength(ncdims);j++) {
	    CDFnode* dim = (CDFnode*)nclistget(ncdims,j);
	    if(dim->dim.declsize == 0) {
	 	/* make node invisible */
		var->invisible = 1;
		var->zerodim = 1;
	    }
	}
    }
    return NC_NOERR;
}

static void
applyclientparamcontrols(NCDAPCOMMON* dapcomm)
{
    /* clear the flags */
    CLRFLAG(dapcomm->controls,NCF_CACHE);
    CLRFLAG(dapcomm->controls,NCF_SHOWFETCH);
    CLRFLAG(dapcomm->controls,NCF_NC3);
    CLRFLAG(dapcomm->controls,NCF_NCDAP);
    CLRFLAG(dapcomm->controls,NCF_PREFETCH);
    CLRFLAG(dapcomm->controls,NCF_PREFETCH_EAGER);

    /* Turn on any default on flags */
    SETFLAG(dapcomm->controls,DFALT_ON_FLAGS);    
    SETFLAG(dapcomm->controls,(NCF_NC3|NCF_NCDAP));

    /* enable/disable caching */
    if(dapparamcheck(dapcomm,"cache",NULL))
	SETFLAG(dapcomm->controls,NCF_CACHE);
    else if(dapparamcheck(dapcomm,"nocache",NULL))
	CLRFLAG(dapcomm->controls,NCF_CACHE);

    /* enable/disable cache prefetch and lazy vs eager*/
    if(dapparamcheck(dapcomm,"prefetch","eager")) {
        SETFLAG(dapcomm->controls,NCF_PREFETCH);
        SETFLAG(dapcomm->controls,NCF_PREFETCH_EAGER);
    } else if(dapparamcheck(dapcomm,"prefetch","lazy")
              || dapparamcheck(dapcomm,"prefetch",NULL)) {
        SETFLAG(dapcomm->controls,NCF_PREFETCH);
        CLRFLAG(dapcomm->controls,NCF_PREFETCH_EAGER);
    } else if(dapparamcheck(dapcomm,"noprefetch",NULL))
        CLRFLAG(dapcomm->controls,NCF_PREFETCH);

    if(FLAGSET(dapcomm->controls,NCF_UNCONSTRAINABLE))
	SETFLAG(dapcomm->controls,NCF_CACHE);

    if(dapparamcheck(dapcomm,"show","fetch"))
	SETFLAG(dapcomm->controls,NCF_SHOWFETCH);

    nclog(NCLOGNOTE,"Caching=%d",FLAGSET(dapcomm->controls,NCF_CACHE));

}
