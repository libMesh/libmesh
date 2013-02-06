/*********************************************************************
 *   Copyright 1993, UCAR/Unidata
 *   See netcdf/COPYRIGHT file for copying and redistribution conditions.
 *   $Header: /upc/share/CVS/netcdf-3/libncdap3/getvara3.c,v 1.44 2010/05/27 21:34:08 dmh Exp $
 *********************************************************************/

#include "ncdap3.h"
#include "dapodom.h"
#include "dapdump.h"
#include "ncd3dispatch.h"


#define NEWVARM

/* Define a tracker for memory to support*/
/* the concatenation*/

struct NCMEMORY {
    void* memory;
    char* next; /* where to store the next chunk of data*/
}; 

/* Forward:*/
static NCerror moveto(NCDAPCOMMON*, Getvara*, CDFnode* dataroot, void* memory);
static NCerror movetor(NCDAPCOMMON*, OCdata currentcontent,
		   NClist* path, int depth,
		   Getvara*, int dimindex,
		   struct NCMEMORY*, NClist* segments);

static int findfield(CDFnode* node, CDFnode* subnode);
static int wholeslicepoint(Dapodometer* odom);
static NCerror removepseudodims(DCEprojection* proj);

static int extract(NCDAPCOMMON*, Getvara*, CDFnode*, DCEsegment*, OClink, OCdata, struct NCMEMORY*);
static int extractstring(NCDAPCOMMON*, Getvara*, CDFnode*, DCEsegment*, OClink, OCdata, struct NCMEMORY*);

/**
1. We build the projection to be sent to the server aka
the fetch constraint.  We want the server do do as much work
as possible, so we send it a url with a fetch constraint
that is the merge of the url constraint with the vara
constraint.

The url constraint, if any, is the one that was provided
in the url specified in nc_open().

The vara constraint is the one formed from the arguments
(start, count, stride) provided to the call to nc_get_vara().

There are some exceptions to the formation of the fetch constraint.
In all cases, the fetch constraint will use any URL selections,
but will use different fetch projections.
a. URL is unconstrainable (e.g. file://...):
	   fetchprojection = null => fetch whole dataset
b. The target variable (as specified in nc_get_vara())
   is already in the cache and is whole variable.
	   fetchprojection = N.A. since variable is in the cache
c. Vara is requesting part of a variable but NCF_WHOLEVAR flag is set.
	   fetchprojection = unsliced vara variable => fetch whole variable
d. Vara is requesting part of a variable and NCF_WHOLEVAR flag is not set.
	   fetchprojection = sliced vara variable => fetch part variable

2. At this point, all or part of the target variable is available in the cache.

3. We build a projection to walk (guide) the use of the oc
   data procedures in extract the required data from the cache.
   For cases a,b,c:
       walkprojection = merge(urlprojection,varaprojection)
   For case d:
       walkprojection =  varaprojection without slicing.
       This means we need only extract the complete contents of the cache.
       Notice that this will not necessarily be a direct memory to
       memory copy because the dap encoding still needs to be
       interpreted. For this case, we derive a walk projection
       from the vara projection that will properly access the cached data.
       This walk projection shifts the merged projection so all slices
       start at 0 and have a stride of 1.

*/

NCerror
nc3d_getvarx(int ncid, int varid,
	    const size_t *startp,
	    const size_t *countp,
	    const ptrdiff_t* stridep,
	    void *data,
	    nc_type dsttype0)
{
    NCerror ncstat = NC_NOERR;
    OCerror ocstat = OC_NOERR;
    int i;
    NC* drno;
    NC* substrate;
    NCDAPCOMMON* dapcomm;
    CDFnode* cdfvar; /* cdf node mapping to var*/
    NClist* varnodes;
    nc_type dsttype;
    Getvara* varainfo = NULL;
    CDFnode* xtarget = NULL; /* target in DATADDS */
    CDFnode* target = NULL; /* target in constrained DDS */
    DCEprojection* varaprojection = NULL;
    NCcachenode* cachenode = NULL;
    size_t localcount[NC_MAX_VAR_DIMS];
    NClist* ncdimsall;
    size_t ncrank;
    NClist* vars = NULL;
    DCEconstraint* fetchconstraint = NULL;
    DCEprojection* fetchprojection = NULL;
    DCEprojection* walkprojection = NULL;
    int state;
#define FETCHWHOLE 1 /* fetch whole data set */
#define FETCHVAR   2 /* fetch whole variable */
#define FETCHPART  4 /* fetch constrained variable */
#define CACHED     8 /* whole variable is already in the cache */

    ncstat = NC_check_id(ncid, (NC**)&drno); 
    if(ncstat != NC_NOERR) goto fail;
    dapcomm = (NCDAPCOMMON*)drno->dispatchdata;
    
    ncstat = NC_check_id(drno->substrate, (NC**)&substrate); 
    if(ncstat != NC_NOERR) goto fail;

    /* Locate var node via varid */
    varnodes = dapcomm->cdf.varnodes;
    for(i=0;i<nclistlength(varnodes);i++) {
	CDFnode* node = (CDFnode*)nclistget(varnodes,i);
	if(node->array.basevar == NULL
           && node->nctype == NC_Primitive
           && node->ncid == varid) {
	    cdfvar = node;
	    break;
	}
    }

    ASSERT((cdfvar != NULL));

    /* Get the dimension info */
    ncdimsall = cdfvar->array.dimsetall;
    ncrank = nclistlength(ncdimsall);

#ifdef DEBUG
 {
int i;
fprintf(stderr,"getvarx: %s",cdfvar->ncfullname);
for(i=0;i<ncrank;i++)
  fprintf(stderr,"[%ld:%ld:%ld]",
	(long)startp[i],
	(long)countp[i],
	(long)stridep[i]
	);
fprintf(stderr,"\n");
 }
#endif

    /* Fill in missing arguments */
    if(startp == NULL)
	startp = nc_sizevector0;

    if(countp == NULL) {
        /* Accumulate the dimension sizes */
        for(i=0;i<ncrank;i++) {
	    CDFnode* dim = (CDFnode*)nclistget(ncdimsall,i);
	    localcount[i] = dim->dim.declsize;
	}
	countp = localcount;
    }

    if(stridep == NULL)
	stridep = nc_ptrdiffvector1;

    /* Validate the dimension sizes */
    for(i=0;i<ncrank;i++) {
        CDFnode* dim = (CDFnode*)nclistget(ncdimsall,i);
	if(startp[i] > dim->dim.declsize
	   || startp[i]+countp[i] > dim->dim.declsize) {
	    ncstat = NC_EINVALCOORDS;
	    goto fail;	    
	}
    }	     

#ifdef DEBUG
 {
NClist* dims = cdfvar->array.dimsetall;
fprintf(stderr,"getvarx: %s",cdfvar->ncfullname);
if(nclistlength(dims) > 0) {int i;
for(i=0;i<nclistlength(dims);i++) 
fprintf(stderr,"[%lu:%lu:%lu]",(unsigned long)startp[i],(unsigned long)countp[i],(unsigned long)stridep[i]);
fprintf(stderr," -> ");
for(i=0;i<nclistlength(dims);i++) 
if(stridep[i]==1)
fprintf(stderr,"[%lu:%lu]",(unsigned long)startp[i],(unsigned long)((startp[i]+countp[i])-1));
else
fprintf(stderr,"[%lu:%lu:%lu]",
(unsigned long)startp[i],
(unsigned long)stridep[i],
(unsigned long)(((startp[i]+countp[i])*stridep[i])-1));
}
fprintf(stderr,"\n");
 }
#endif

    dsttype = (dsttype0);

    /* Default to using the inquiry type for this var*/
    if(dsttype == NC_NAT) dsttype = cdfvar->externaltype;

    /* Validate any implied type conversion*/
    if(cdfvar->etype != dsttype && dsttype == NC_CHAR) {
        /* The only disallowed conversion is to/from char and non-byte
           numeric types*/
	switch (cdfvar->etype) {
	case NC_STRING: case NC_URL:
	case NC_CHAR: case NC_BYTE: case NC_UBYTE:
	    break;
	default:
	    return THROW(NC_ECHAR);
	}
    }

    ncstat = makegetvar34(dapcomm,cdfvar,data,dsttype,&varainfo);
    if(ncstat != NC_NOERR) {THROWCHK(ncstat); goto fail;}

    state = 0;
    if(iscached(dapcomm,cdfvar,&cachenode)) {
	state = CACHED;
	ASSERT((cachenode != NULL));
#ifdef DEBUG
fprintf(stderr,"var is in cache\n");
#endif
        /* If it is cached, then it is a whole variable but may still
           need to apply constraints during the walk */
	ASSERT(cachenode->wholevariable); /* by construction */
    } else if(FLAGSET(dapcomm->controls,NCF_UNCONSTRAINABLE)) {
	state = FETCHWHOLE;
    } else {/* load using constraints */
        if(FLAGSET(dapcomm->controls,NCF_WHOLEVAR))
	    state = FETCHVAR;
	else
	    state = FETCHPART;
    }

    ASSERT(state != 0);    

    /* Convert the start/stop/stride info into a projection */
    ncstat = buildvaraprojection3(varainfo,
		                  startp,countp,stridep,
                                  &varaprojection);
    if(ncstat != NC_NOERR) {THROWCHK(ncstat); goto fail;}

    fetchprojection = NULL;
    walkprojection = NULL;

    /* Create walkprojection as the merge of the url projections
       and the vara projection; may change in FETCHPART case below*/
    ncstat = daprestrictprojection(dapcomm->oc.dapconstraint->projections,
				   varaprojection,&walkprojection);
    if(ncstat != NC_NOERR) {THROWCHK(ncstat); goto fail;}

#ifdef DEBUG
fprintf(stderr,"getvarx: walkprojection: |%s|\n",dumpprojection(walkprojection));
#endif

    /* define the var list of interest */
    vars = nclistnew();
    nclistpush(vars,(ncelem)varainfo->target);

    switch (state) {

    case FETCHWHOLE: {
        /* buildcachenode3 will create a new cachenode and
           will also fetch the whole corresponding datadds.
	*/
        /* Build the complete constraint to use in the fetch */
        fetchconstraint = (DCEconstraint*)dcecreate(CES_CONSTRAINT);
        /* Use no projections or selections */
        fetchconstraint->projections = nclistnew();
        fetchconstraint->selections = nclistnew();
#ifdef DEBUG
fprintf(stderr,"getvarx: FETCHWHOLE: fetchconstraint: %s\n",dumpconstraint(fetchconstraint));
#endif
        ncstat = buildcachenode34(dapcomm,fetchconstraint,vars,&cachenode,0);
	fetchconstraint = NULL; /*buildcachenode34 takes control of fetchconstraint.*/
	if(ncstat != NC_NOERR) {THROWCHK(ncstat); goto fail;}
    } break;

    case CACHED: {
    } break;

    case FETCHVAR: { /* Fetch a complete single variable */
        /* Create fetch projection as the merge of the url projections
           and the vara projection */
        ncstat = daprestrictprojection(dapcomm->oc.dapconstraint->projections,
				       varaprojection,&fetchprojection);
	/* elide any sequence and string dimensions (dap servers do not allow such). */
	ncstat = removepseudodims(fetchprojection);
        if(ncstat != NC_NOERR) {THROWCHK(ncstat); goto fail;}

	/* Convert to a whole variable projection */
	dcemakewholeprojection(fetchprojection);

#ifdef DEBUG
fprintf(stderr,"getvarx: FETCHVAR: fetchprojection: |%s|\n",dumpprojection(fetchprojection));
#endif

        /* Build the complete constraint to use in the fetch */
        fetchconstraint = (DCEconstraint*)dcecreate(CES_CONSTRAINT);
        /* merged constraint just uses the url constraint selection */
        fetchconstraint->selections = dceclonelist(dapcomm->oc.dapconstraint->selections);
	/* and the created fetch projection */
        fetchconstraint->projections = nclistnew();
	nclistpush(fetchconstraint->projections,(ncelem)fetchprojection);
#ifdef DEBUG
fprintf(stderr,"getvarx: FETCHVAR: fetchconstraint: %s\n",dumpconstraint(fetchconstraint));
#endif
        /* buildcachenode3 will create a new cachenode and
           will also fetch the corresponding datadds.
        */
        ncstat = buildcachenode34(dapcomm,fetchconstraint,vars,&cachenode,0);
	fetchconstraint = NULL; /*buildcachenode34 takes control of fetchconstraint.*/
	if(ncstat != NC_NOERR) {THROWCHK(ncstat); goto fail;}
    } break;

    case FETCHPART: {
        /* Create fetch projection as the merge of the url projections
           and the vara projection */
        ncstat = daprestrictprojection(dapcomm->oc.dapconstraint->projections,
				       varaprojection,&fetchprojection);
	/* elide any sequence and string dimensions (dap servers do not allow such). */
	ncstat = removepseudodims(fetchprojection);
        if(ncstat != NC_NOERR) {THROWCHK(ncstat); goto fail;}

	/* Shift the varaprojection for simple walk */
	dcefree((DCEnode*)walkprojection) ; /* reclaim any existing walkprojection */        
	walkprojection = (DCEprojection*)dceclone((DCEnode*)varaprojection);
        dapshiftprojection(walkprojection);

#ifdef DEBUG
fprintf(stderr,"getvarx: FETCHPART: fetchprojection: |%s|\n",dumpprojection(fetchprojection));
#endif

        /* Build the complete constraint to use in the fetch */
        fetchconstraint = (DCEconstraint*)dcecreate(CES_CONSTRAINT);
        /* merged constraint just uses the url constraint selection */
        fetchconstraint->selections = dceclonelist(dapcomm->oc.dapconstraint->selections);
	/* and the created fetch projection */
        fetchconstraint->projections = nclistnew();
	nclistpush(fetchconstraint->projections,(ncelem)fetchprojection);
#ifdef DEBUG
fprintf(stderr,"getvarx: FETCHPART: fetchconstraint: %s\n",dumpconstraint(fetchconstraint));
#endif
        /* buildcachenode3 will create a new cachenode and
           will also fetch the corresponding datadds.
        */
        ncstat = buildcachenode34(dapcomm,fetchconstraint,vars,&cachenode,0);
	fetchconstraint = NULL; /*buildcachenode34 takes control of fetchconstraint.*/
	if(ncstat != NC_NOERR) {THROWCHK(ncstat); goto fail;}
    } break;

    default: PANIC1("unknown fetch state: %d\n",state);
    }

    ASSERT(cachenode != NULL);

#ifdef DEBUG
fprintf(stderr,"cache.datadds=%s\n",dumptree(cachenode->datadds));
#endif

    /* attach DATADDS to (constrained) DDS */
    unattach34(dapcomm->cdf.ddsroot);
    ncstat = attachsubset34(cachenode->datadds,dapcomm->cdf.ddsroot);
    if(ncstat) goto fail;	

    /* Fix up varainfo to use the cache */
    varainfo->cache = cachenode;
    cachenode = NULL;
    varainfo->varaprojection = walkprojection;
    walkprojection = NULL;

    /* Get the var correlate from the datadds */
    target = varainfo->target;
    xtarget = target->attachment;
    if(xtarget == NULL) 
	{THROWCHK(ncstat=NC_ENODATA); goto fail;}

    /* Switch to datadds tree space*/
    varainfo->target = xtarget;
    ncstat = moveto(dapcomm,varainfo,varainfo->cache->datadds,data);
    if(ncstat != NC_NOERR) {THROWCHK(ncstat); goto fail;}

    goto ok;

fail:
    if(ocstat != OC_NOERR) ncstat = ocerrtoncerr(ocstat);
ok:
    nclistfree(vars);
    dcefree((DCEnode*)varaprojection);
    dcefree((DCEnode*)fetchconstraint);
    freegetvara(varainfo);
    return THROW(ncstat);
}

/* Remove any pseudodimensions (sequence and string)*/
static NCerror
removepseudodims(DCEprojection* proj)
{
    int i;
#ifdef DEBUG1
fprintf(stderr,"removesequencedims.before: %s\n",dumpprojection(proj));
#endif
    for(i=0;i<nclistlength(proj->var->segments);i++) {
	DCEsegment* seg = (DCEsegment*)nclistget(proj->var->segments,i);
	CDFnode* cdfnode = (CDFnode*)seg->annotation;
	if(cdfnode->array.seqdim != NULL)
	    seg->rank = 0;
	else if(cdfnode->array.stringdim != NULL)
	    seg->rank--;
    }
#ifdef DEBUG1
fprintf(stderr,"removepseudodims.after: %s\n",dumpprojection(proj));
#endif
    return NC_NOERR;
}

static NCerror
moveto(NCDAPCOMMON* nccomm, Getvara* xgetvar, CDFnode* xrootnode, void* memory)
{
    OCerror ocstat = OC_NOERR;
    NCerror ncstat = NC_NOERR;
    OCconnection conn = nccomm->oc.conn;
    OCdata xrootcontent;
    OCobject ocroot;
    NClist* path = nclistnew();
    struct NCMEMORY memstate;

    memstate.next = (memstate.memory = memory);

    /* Get the root content*/
    ocroot = xrootnode->tree->ocroot;
    xrootcontent = oc_data_new(conn);
    ocstat = oc_data_root(conn,ocroot,xrootcontent);
    if(ocstat != OC_NOERR) {THROWCHK(ocstat); goto done;}

    /* Remember: xgetvar->target is in DATADDS tree */
    collectnodepath3(xgetvar->target,path,WITHDATASET);
    ncstat = movetor(nccomm,xrootcontent,
                     path,0,xgetvar,0,&memstate,
                     xgetvar->varaprojection->var->segments);

done:
    nclistfree(path);
    oc_data_free(conn,xrootcontent);
    if(ocstat != OC_NOERR) ncstat = ocerrtoncerr(ocstat);
    return THROW(ncstat);
}

static NCerror
movetor(NCDAPCOMMON* nccomm,
	OCdata currentcontent,
	NClist* path,
        int depth, /* depth is position in segment list*/
	Getvara* xgetvar,
        int dimindex, /* dimindex is position in xgetvar->slices*/
	struct NCMEMORY* memory,
	NClist* segments)
{
    int i;
    OCerror ocstat = OC_NOERR;
    NCerror ncstat = NC_NOERR;
    size_t fieldindex,gridindex,rank;
    OCconnection conn = nccomm->oc.conn;
    CDFnode* xnode = (CDFnode*)nclistget(path,depth);
    OCdata reccontent = OCNULL;
    OCdata dimcontent = OCNULL;
    OCdata fieldcontent = OCNULL;
    Dapodometer* odom = OCNULL;
    OCmode currentmode = OCNULLMODE;
    CDFnode* xnext;
    int hasstringdim = 0;
    size_t dimoffset;
    DCEsegment* segment;
    int newdepth;
    int caching = FLAGSET(nccomm->controls,NCF_CACHE);
    int unconstrainable = FLAGSET(nccomm->controls,NCF_UNCONSTRAINABLE);

    /* Note that we use depth-1 because the path contains the DATASET
       but the segment list does not */
    segment = (DCEsegment*)nclistget(segments,depth-1); /*may be NULL*/
    if(xnode->etype == NC_STRING || xnode->etype == NC_URL) hasstringdim = 1;

    ocstat = oc_data_mode(conn,currentcontent,&currentmode);

#ifdef DEBUG2
fprintf(stderr,"moveto: nctype=%d currentmode=%d depth=%d dimindex=%d",
        xnode->nctype, currentmode, depth,dimindex);
fprintf(stderr," segment=%s hasstringdim=%d\n",
		dcetostring((DCEnode*)segment),hasstringdim);
#endif

    /* Switch on the combination of nctype and mode */
#define CASE(nc1,nc2) (nc1*1024+nc2)

    /* This must be consistent with the oc mode transition function */
    switch (CASE(xnode->nctype,currentmode)) {

    default:
	PANIC2("Illegal combination: nctype=%d mode=%d",
		(int)xnode->nctype,(int)currentmode);
	break;

    case CASE(NC_Sequence,OCFIELDMODE):
    case CASE(NC_Dataset,OCFIELDMODE):
    case CASE(NC_Grid,OCFIELDMODE):
    case CASE(NC_Structure,OCFIELDMODE):
	/* currentcontent points to the grid/dataset/structure instance */
	xnext = (CDFnode*)nclistget(path,depth+1);
	ASSERT((xnext != NULL));
	fieldindex = findfield(xnode,xnext);
	/* If the next node is a virtual node, then
	   we need to effectively
	   ignore it and use the appropriate subnode.
	   If the next node is a structuregrid node, then
	   use it as is.
	*/
        if(xnext->virtual) {
	    CDFnode* xgrid = xnext;
	    xnext = (CDFnode*)nclistget(path,depth+2); /* real node */
	    gridindex = fieldindex;
	    fieldindex = findfield(xgrid,xnext);
	    fieldindex += gridindex;
	    newdepth = depth+2;
	} else {
	    newdepth = depth+1;
	}
        fieldcontent = oc_data_new(conn);
        ocstat = oc_data_ith(conn,currentcontent,fieldindex,fieldcontent);
	if(ocstat != OC_NOERR) {THROWCHK(ocstat); goto fail;}
	ncstat = movetor(nccomm,fieldcontent,
                         path,newdepth,xgetvar,dimindex,memory,
			 segments);
	break;

    case CASE(NC_Sequence,OCARRAYMODE): /* will actually always be scalar, but will have
                                           rank == 1 to account for the sequence dim */
    case CASE(NC_Grid,OCARRAYMODE): /* will actually always be scalar */
    case CASE(NC_Structure,OCARRAYMODE):
        /* figure out which slices refer to this node:
           dimindex upto dimindex+rank; */
        ASSERT((segment != NULL));
        rank = segment->rank;
	if(xnode->nctype == NC_Sequence)
	    rank--; /* ignore the sequence dim */
	if(rank == 0) {
            odom = newdapodometer1(1);
	} else if(caching || unconstrainable) {	
            odom = newdapodometer(segment->slices,0,rank);	    
	} else { /*Since vara was projected out, build a simple odometer*/
            odom = newsimpledapodometer(segment,rank);
	}
        while(dapodometermore(odom)) {
            OCmode mode;
            /* Compute which instance to move to*/
            dimoffset = dapodometercount(odom);
            dimcontent = oc_data_new(conn);
            ocstat = oc_data_ith(conn,currentcontent,dimoffset,dimcontent);
            if(ocstat != OC_NOERR) {THROWCHK(ocstat); goto fail;}
            ocstat = oc_data_mode(conn,dimcontent,&mode);
            ASSERT((mode == OCFIELDMODE
		    || (mode == OCSEQUENCEMODE && xnode->nctype == NC_Sequence)));
            ncstat = movetor(nccomm,dimcontent,
                                 path,depth,
                                 xgetvar,dimindex+rank,
                                 memory,segments);
            dapodometerincr(odom);
        }
        freedapodometer(odom);
        break;

    case CASE(NC_Sequence,OCSEQUENCEMODE): {
        DCEslice* uslice;
        ASSERT((segment != NULL));
        /* Get and check the corresponding sequence dimension from DDS */
        ASSERT((xnode->attachment != NULL));
        /* use uslice to walk the sequence; however, watch out
           for the case when the user set a limit and that limit
           is not actually reached in this request.
        */
        /* By construction, this sequence represents the first 
           (and only) dimension of this segment */
        uslice = &segment->slices[0];
        reccontent = oc_data_new(conn);
        for(i=uslice->first;i<uslice->stop;i+=uslice->stride) {
	    OCmode eos;
            ocstat = oc_data_ith(conn,currentcontent,i,reccontent);
	    if(ocstat != OC_NOERR) {THROWCHK(ocstat); goto fail;}
	    ocstat = oc_data_mode(conn,reccontent,&eos);
	    if(ocstat != OC_NOERR) {THROWCHK(ocstat); goto fail;}
	    if(eos == OCNULLMODE) {
                /* We asked for too much */
                ncstat = THROW(NC_EINVALCOORDS);
                goto fail;
            }
            ncstat = movetor(nccomm,reccontent,
                                 path,depth,
                                 xgetvar,dimindex+1,
                                 memory,segments);
            if(ncstat != OC_NOERR) {THROWCHK(ncstat); goto fail;}
        }
        } break;

    case CASE(NC_Primitive,OCPRIMITIVEMODE):
        if(hasstringdim)
	    ncstat = extractstring(nccomm, xgetvar, xnode, segment, conn, currentcontent, memory);
	else	
	    ncstat = extract(nccomm, xgetvar, xnode, segment, conn, currentcontent, memory);
	break;

    }
    goto ok;

fail:
ok:
    oc_data_free(conn,dimcontent);
    oc_data_free(conn,fieldcontent);
    oc_data_free(conn,reccontent);
    if(ocstat != OC_NOERR) ncstat = ocerrtoncerr(ocstat);
    return THROW(ncstat);
}

/* Determine the index in the odometer at which
   the odometer will be walking the whole subslice
   This will allow us to optimize.
*/
static int 
wholeslicepoint(Dapodometer* odom)
{
    unsigned int i;
    int point;
    for(point=-1,i=0;i<odom->rank;i++) {
        ASSERT((odom->slices[i].declsize != 0));
        if(odom->slices[i].first != 0 || odom->slices[i].stride != 1
           || odom->slices[i].length != odom->slices[i].declsize)
	    point = i;
    }
    if(point == -1)
	point = 0; /* wholevariable */
    else if(point == (odom->rank - 1)) 
	point = -1; /* no whole point */
    else
	point += 1; /* intermediate point */
    return point;
}

static int
findfield(CDFnode* node, CDFnode* field)
{
    size_t i;
    for(i=0;i<nclistlength(node->subnodes);i++) {
        CDFnode* test = (CDFnode*) nclistget(node->subnodes,i);
        if(test == field) return i;
    }
    return -1;
}


int
nc3d_getvarmx(int ncid, int varid,
	    const size_t *start,
	    const size_t *edges,
	    const ptrdiff_t* stride,
 	    const ptrdiff_t* map,
	    void* data,
	    nc_type dsttype0)
{
    NCerror ncstat = NC_NOERR;
    int i;
    NC* drno;
    NC* substrate;
    NCDAPCOMMON* dapcomm;
    NC_var* var;
    CDFnode* cdfvar; /* cdf node mapping to var*/
    NClist* varnodes;
    nc_type dsttype;
    size_t externsize;
    size_t dimsizes[NC_MAX_VAR_DIMS];
    Dapodometer* odom = NULL;
    unsigned int ncrank;
    NClist* ncdims = NULL;
    size_t nelems;
#ifdef NEWVARM
    char* localcopy; /* of whole variable */
#endif

    ncstat = NC_check_id(ncid, (NC**)&drno); 
    if(ncstat != NC_NOERR) goto done;
    dapcomm = (NCDAPCOMMON*)drno->dispatchdata;

    ncstat = NC_check_id(drno->substrate, (NC**)&substrate); 
    if(ncstat != NC_NOERR) goto done;
    var = NC_lookupvar(substrate,varid);
    if(var == NULL) {ncstat = NC_ENOTVAR; goto done;}

    /* Locate var node via varid */
    varnodes = dapcomm->cdf.varnodes;
    for(i=0;i<nclistlength(varnodes);i++) {
	CDFnode* node = (CDFnode*)nclistget(varnodes,i);
	if(node->array.basevar == NULL
           && node->nctype == NC_Primitive
           && node->ncid == varid) {
	    cdfvar = node;
	    break;
	}
    }

    ASSERT((cdfvar != NULL));
    ASSERT((strcmp(cdfvar->ncfullname,var->name->cp)==0));

    if(nclistlength(cdfvar->array.dimsetplus) == 0) {
       /* The variable is a scalar; consequently, there is only one
          thing to get and only one place to put it.  (Why was I
          called?) */
	/* recurse with additional parameters */
        return THROW(nc3d_getvarx(ncid,varid,
		 NULL,NULL,NULL,
		 data,dsttype0));
    }
         
    dsttype = (dsttype0);

    /* Default to using the inquiry type for this var*/
    if(dsttype == NC_NAT) dsttype = cdfvar->externaltype;

    /* Validate any implied type conversion*/
    if(cdfvar->etype != dsttype && dsttype == NC_CHAR) {
	/* The only disallowed conversion is to/from char and non-byte
           numeric types*/
	switch (cdfvar->etype) {
	case NC_STRING: case NC_URL:
	case NC_CHAR: case NC_BYTE: case NC_UBYTE:
 	    break;
	default:
	    return THROW(NC_ECHAR);
	}
    }

    externsize = nctypesizeof(dsttype);

    /* Accumulate the dimension sizes and the total # of elements */
    ncdims = cdfvar->array.dimsetall;
    ncrank = nclistlength(ncdims);

    nelems = 1; /* also Compute the number of elements being retrieved */
    for(i=0;i<ncrank;i++) {
	CDFnode* dim = (CDFnode*)nclistget(ncdims,i);
	dimsizes[i] = dim->dim.declsize;
	nelems *= edges[i];
    }

    /* Originally, this code repeatedly extracted single values
       using get_var1. In an attempt to improve performance,
       I have converted to reading the whole variable at once
       and walking it locally.
    */

#ifdef NEWVARM
    localcopy = (char*)malloc(nelems*externsize);

    /* We need to use the varieties of get_vars in order to
       properly do conversion to the external type
    */

    switch (dsttype) {

    case NC_CHAR:
	ncstat = nc_get_vars_text(ncid,varid,start, edges, stride,
				  (char*)localcopy);
	break;
    case NC_BYTE:
	ncstat = nc_get_vars_schar(ncid,varid,start, edges, stride,
				   (signed char*)localcopy);
	break;
    case NC_SHORT:
	ncstat = nc_get_vars_short(ncid,varid, start, edges, stride,
			  	   (short*)localcopy);
	break;
    case NC_INT:
	ncstat = nc_get_vars_int(ncid,varid,start, edges, stride,
				 (int*)localcopy);
	break;
    case NC_FLOAT:
	ncstat = nc_get_vars_float(ncid,varid,start, edges, stride,
				   (float*)localcopy);
	break;
    case NC_DOUBLE:
	ncstat = nc_get_vars_double(ncid,varid,	start, edges, stride,
		 		    (double*)localcopy);
	break;
    default: break;
    }

    odom = newdapodometer2(start,edges,stride,0,ncrank);

    /* Walk the local copy */
    for(i=0;i<nelems;i++) {
	size_t voffset = dapodometervarmcount(odom,map,dimsizes);
	void* dataoffset = (void*)(((char*)data) + (externsize*voffset));
	char* localpos = (localcopy + externsize*i);
	/* extract the indexset'th value from local copy */
	memcpy(dataoffset,(void*)localpos,externsize);
/*
fprintf(stderr,"new: %lu -> %lu  %f\n",
	(unsigned long)(i),
        (unsigned long)voffset,
	*(float*)localpos);
*/
	dapodometerincr(odom);
    }    
#else
    odom = newdapodometer2(start,edges,stride,0,ncrank);
    while(dapodometermore(odom)) {
	size_t* indexset = dapodometerindices(odom);
	size_t voffset = dapodometervarmcount(odom,map,dimsizes);
	char internalmem[128];
	char externalmem[128];
	void* dataoffset = (void*)(((char*)data) + (externsize*voffset));

	/* get the indexset'th value using variable's internal type */
	ncstat = nc_get_var1(ncid,varid,indexset,(void*)&internalmem);
        if(ncstat != NC_NOERR) goto done;
	/* Convert to external type */
	ncstat = dapconvert3(cdfvar->etype,dsttype,externalmem,internalmem);
        if(ncstat != NC_NOERR) goto done;
	memcpy(dataoffset,(void*)externalmem,externsize);
/*
fprintf(stderr,"old: %lu -> %lu  %f\n",
	(unsigned long)dapodometercount(odom),
        (unsigned long)voffset,
	*(float*)externalmem);
*/
	dapodometerincr(odom);
    }    
#endif

done:
    return ncstat;
}

static int
conversionrequired(nc_type t1, nc_type t2)
{
    if(t1 == t2)
	return 0;
    if(nctypesizeof(t1) != nctypesizeof(t2))
	return 1;
    /* Avoid too many cases by making t1 < t2 */
    if(t1 > t2) {int tmp = t1; t1 = t2; t2 = tmp;}
#undef CASE
#define CASE(t1,t2) ((t1)<<5 | (t2))
    switch (CASE(t1,t2)) {
    case CASE(NC_BYTE,NC_UBYTE):
    case CASE(NC_BYTE,NC_CHAR):
    case CASE(NC_CHAR,NC_UBYTE):
    case CASE(NC_SHORT,NC_USHORT):
    case CASE(NC_INT,NC_UINT):
    case CASE(NC_INT64,NC_UINT64):
	return 0;
    default: break;
    }
    return 1;
}

/* We are at a primitive variable or scalar that has no string dimensions.
Extract the data.
(This is way too complicated)
*/
static int
extract(
	NCDAPCOMMON* nccomm,
	Getvara* xgetvar,
	CDFnode* xnode,
        DCEsegment* segment,
        OClink conn,
        OCdata currentcontent,
	struct NCMEMORY* memory
       )
{
    OCerror ocstat = OC_NOERR;
    NCerror ncstat = NC_NOERR;
    size_t rank;
    Dapodometer* odom = OCNULL;
    int wholepoint;
    size_t externtypesize;
    size_t interntypesize;
    char* localmemory = NULL;
    size_t odomsubsize;
    size_t internlen;
    int requireconversion;
    char value[16]; 

    ASSERT((segment != NULL));

    requireconversion = conversionrequired(xgetvar->dsttype,xnode->etype);

    rank = segment->rank;

    if(rank == 0) {/* scalar */
	char* mem = (requireconversion?value:memory->next);
        ASSERT((segment != NULL));
        externtypesize = nctypesizeof(xgetvar->dsttype);
	ASSERT(externtypesize <= sizeof(value));
	/* Read the whole scalar directly into memory  */
	ocstat = oc_data_get(conn,currentcontent,mem,externtypesize,0,1);
	if(ocstat != OC_NOERR) {THROWCHK(ocstat); goto done;}
	if(requireconversion) {
	    /* convert the value to external type */
            ncstat = dapconvert3(xnode->etype,xgetvar->dsttype,memory->next,value,1);
            if(ncstat != NC_NOERR) {THROWCHK(ncstat); goto done;}
        }
        memory->next += (externtypesize);

    } else {/* rank > 0 */

#ifdef DEBUG2
fprintf(stderr,"moveto: primitive: segment=%s",
                dcetostring((DCEnode*)segment));
fprintf(stderr," iswholevariable=%d",xgetvar->cache->wholevariable);
fprintf(stderr,"\n");
#endif

        ASSERT(xgetvar->cache != NULL);
        if(xgetvar->cache->wholevariable) {
            odom = newdapodometer(segment->slices,0,rank);
        } else { /*!xgetvar->cache->wholevariable*/
            odom = newsimpledapodometer(segment,rank);
        }
        /* Optimize off the use of the odometer by checking the slicing
           to see if the whole variable, or some whole subslice
           is being extracted.
           However do not do this if the external type conversion is needed
           or if the whole slice point is rank-1 (normal case anyway).
        */
        externtypesize = nctypesizeof(xgetvar->dsttype);
        interntypesize = nctypesizeof(xnode->etype);
        wholepoint = wholeslicepoint(odom);
        if(wholepoint == -1)
            odomsubsize = 1; /* no whole point */
        else
            odomsubsize = dapodometerspace(odom,wholepoint);
        internlen = (odomsubsize*interntypesize);
        if(requireconversion) {
            /* copy the data locally before conversion */
            localmemory = (char*)malloc(internlen);
        } else {
            localmemory = memory->next;
        }

#ifdef DEBUG2
fprintf(stderr,"moveto: primitive: ");
fprintf(stderr," wholepoint=%d",wholepoint);
fprintf(stderr,"\n");
#endif

        if(wholepoint == 0) {/* whole variable */
            /* Read the whole n elements directly into memory.*/
            ocstat = oc_data_get(conn,currentcontent,localmemory,
                                 internlen,0,odomsubsize);
            if(ocstat != OC_NOERR) {THROWCHK(ocstat); goto done;}
	    if(requireconversion) {
                /* do conversion */
                ncstat = dapconvert3(xnode->etype,xgetvar->dsttype,
                                     memory->next,localmemory,odomsubsize);
                if(ncstat != NC_NOERR) {THROWCHK(ncstat); goto done;}
            }
            memory->next += (externtypesize*odomsubsize);
        } else if(wholepoint > 0) {/* whole subslice */
            odom->rank = wholepoint; /* truncate */
            while(dapodometermore(odom)) {
                size_t dimoffset = dapodometercount(odom) * odomsubsize;
                ocstat = oc_data_get(conn,currentcontent,localmemory,
                                         internlen,dimoffset,odomsubsize);
                if(ocstat != OC_NOERR) {THROWCHK(ocstat); goto done;}
		if(requireconversion) {
                    /* do conversion */
                    ncstat = dapconvert3(xnode->etype,xgetvar->dsttype,
                                         memory->next,localmemory,odomsubsize);
                    if(ncstat != NC_NOERR) {THROWCHK(ncstat); goto done;}
                }
                memory->next += (externtypesize*odomsubsize);
                dapodometerincr(odom);
            }
        } else { /* Oh well, use the odometer to walk to the
                    appropriate fields*/
            while(dapodometermore(odom)) {
		char* mem = (requireconversion?value:memory->next);
                size_t dimoffset = dapodometercount(odom);
                ocstat = oc_data_get(conn,currentcontent,mem,externtypesize,dimoffset,1);
                if(ocstat != OC_NOERR) {THROWCHK(ocstat); goto done;}
		if(requireconversion) {
                    ncstat = dapconvert3(xnode->etype,xgetvar->dsttype,memory->next,value,1);
                    if(ncstat != NC_NOERR) {THROWCHK(ncstat); goto done;}
		}
                memory->next += externtypesize;
                dapodometerincr(odom);
            }
        }
        freedapodometer(odom);
        if(requireconversion) nullfree(localmemory);
    }
done:
    return THROW(ncstat);
}


static NCerror
slicestring(OCconnection conn, char* stringmem, DCEslice* slice, struct NCMEMORY* memory)
{
    size_t stringlen;
    unsigned int i;
    NCerror ncstat = NC_NOERR;
    char* lastchar;
    size_t charcount; /* number of characters inserted into memory */

    /* libnc-dap chooses to convert string escapes to the corresponding
       character; so we do likewise.
    */
    dapexpandescapes(stringmem); 
    stringlen = strlen(stringmem);

#ifdef DEBUG2
fprintf(stderr,"moveto: slicestring: string/%lu=%s\n",stringlen,stringmem);
fprintf(stderr,"slicestring: %lu string=|%s|\n",stringlen,stringmem);
fprintf(stderr,"slicestring: slice=[%lu:%lu:%lu/%lu]\n",
slice->first,slice->stride,slice->stop,slice->declsize);
#endif

    /* Stride across string; if we go past end of string, then pad*/
    charcount = 0;
    for(i=slice->first;i<slice->length;i+=slice->stride) {
        if(i < stringlen)
            *memory->next = stringmem[i];
        else /* i >= stringlen*/
            *memory->next = NC_FILL_CHAR;
	memory->next++;
	charcount++;
    }
    lastchar = (memory->next);
    if(charcount > 0) {
        lastchar--;
    }

    return THROW(ncstat);
}

/*
Extract data for a netcdf variable that has a string dimension.
*/
static int
extractstring(
	NCDAPCOMMON* nccomm,
	Getvara* xgetvar,
	CDFnode* xnode,
        DCEsegment* segment,
        OClink conn,
        OCdata currentcontent,
	struct NCMEMORY* memory
       )
{
    NCerror ncstat = NC_NOERR;
    OCerror ocstat = OC_NOERR;
    int i;
    size_t rank;
    int caching = FLAGSET(nccomm->controls,NCF_CACHE);
    int unconstrainable = FLAGSET(nccomm->controls,NCF_UNCONSTRAINABLE);
    NClist* strings = NULL;
    Dapodometer* odom = OCNULL;

    rank = segment->rank;

    /* A number of optimizations are possible but none is currently used. */

    /* Use the odometer to walk to the appropriate fields*/
    if(rank == 1) {
        odom = newdapodometer1(1); /* scalar case */
    } else if(caching || unconstrainable) {	
        odom = newdapodometer(segment->slices,0,rank-1);	    
    } else { /*Since vara was projected out, build a simple odometer*/
        odom = newsimpledapodometer(segment,rank-1);
    }

    /* step thru the odometer obtaining each string and storing it in an OClist */
    strings = nclistnew();
    nclistsetalloc(strings,dapodometerspace(odom,0)); /* preallocate */
    while(dapodometermore(odom)) {
	char* value = NULL;
	size_t dimoffset = dapodometercount(odom);
	ocstat = oc_data_get(conn,currentcontent,&value,sizeof(value),dimoffset,1);
	if(ocstat != OC_NOERR) goto done;
	nclistpush(strings,(ncelem)value);	
        dapodometerincr(odom);
    }
    freedapodometer(odom);
    /* Get each string in turn, slice it and store in user
       supplied memory */
    for(i=0;i<nclistlength(strings);i++) {
	char* s = (char*)nclistget(strings,i);
	slicestring(conn,s,&segment->slices[rank-1],memory);
	free(s);	
    }    
    nclistfree(strings);
done:
    if(ocstat != OC_NOERR) ncstat = ocerrtoncerr(ocstat);
    return THROW(ncstat);
}
