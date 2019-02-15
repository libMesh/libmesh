/*********************************************************************
 *   Copyright 2010, UCAR/Unidata
 *   See netcdf/COPYRIGHT file for copying and redistribution conditions.
 *   $Id $
 *********************************************************************/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "netcdf.h"
#include "list.h"
#include "utils.h"
#include "chunkspec.h"

/* Structure mapping dimension IDs to corresponding chunksizes. */
static struct DimChunkSpecs {
    size_t ndims;		/* number of dimensions in chunkspec string */
    int *idimids;		/* (input) ids for dimensions in chunkspec string */
    size_t *chunksizes;		/* corresponding chunk sizes */
    bool_t omit;		/* true if chunking to be turned off */
} dimchunkspecs;

struct VarChunkSpec {
    size_t rank;			/* number of dimensions in chunkspec string */
    size_t chunksizes[NC_MAX_VAR_DIMS]; /* corresponding chunk sizes */
    bool_t omit;			/* true if chunking to be turned off */
    int igrpid;				/* container of the (input) variable */
    int ivarid;				/* (input) Variable whose chunks are specified */
};

static List* varchunkspecs = NULL;      /* List<VarChunkSpec> */

/* Forward */
static int dimchunkspec_parse(int ncid, const char *spec);
static int varchunkspec_parse(int ncid, const char *spec);

void
chunkspecinit(void)
{
    /* initialization */
    if(varchunkspecs == NULL)
	varchunkspecs = listnew();
    memset(&dimchunkspecs,0,sizeof(dimchunkspecs));
}

/*
 * Parse chunkspec string of either kind.
 * Returns NC_NOERR if no error, NC_EINVAL if spec was malformed.
 */
int
chunkspec_parse(int igrp, const char *spec)
{
    /* Decide if this is a per-variable or per-dimension chunkspec */
    if (!spec || *spec == '\0')
	return NC_NOERR; /* Use defaults */
    if(strchr(spec,':') == NULL)
	return dimchunkspec_parse(igrp,spec);
    else
	return varchunkspec_parse(igrp,spec);
}

/*
 * Parse chunkspec string and convert into dimchunkspec structure.
 *   ncid: location ID of open netCDF file or group in an open file
 *   spec: string of form 
 *           dim1/n1,dim2/n2,...,dimk/nk
 *   specifying chunk size (ni) to be used for dimension named
 *   dimi.  Dimension names may be absolute,
 *   e.g. "/grp_a/grp_a1/dim".  The "ni" part of the spec may be
 *   omitted, in which case it is assumed to be the entire
 *   dimension size.  That is also the default for dimensions
 *   not mentioned in the string. However, for unlimited dimensions,
 *   the default is a default size: 4 megabytes or the
 *   existing unlimited size if smaller.
 *   If the chunkspec string is "/", specifying no dimensions or 
 *   chunk sizes, it indicates chunking to be turned off on output.
 *
 * Returns NC_NOERR if no error, NC_EINVAL if spec has consecutive
 * unescaped commas or no chunksize specified for dimension.
 */
static int
dimchunkspec_parse(int igrp, const char *spec)
{
    const char *cp;	   /* character cursor */
    const char *pp = spec; /* previous char cursor for detecting escapes */
    const char *np;	   /* beginning of current dimension name */
    size_t ndims = 0;
    int idim;
    int ret = NC_NOERR;
    int comma_seen = 0;

    dimchunkspecs.ndims = 0;
    dimchunkspecs.omit = false;
    if (!spec || *spec == '\0') /* default chunking */
	goto done;
    /* Special rule: // is treated as equivalent to / */
    if ((spec[0] == '/' && spec[1] == '\0')
	|| (spec[0] == '/' && spec[1] == '/' && spec[2] == '\0')) { /* no chunking */
	dimchunkspecs.omit = true;
	goto done;
    }
    /* Count unescaped commas, handle consecutive unescaped commas as error */
    for(cp = spec; *cp; cp++) {
	if(*cp == ',' && *pp != '\\') {
	    if(comma_seen) {	/* consecutive commas detected */
		{ret = NC_EINVAL; goto done;}
	    }
	    comma_seen = 1;
	    ndims++;
	} else {
	    comma_seen = 0;
	}
	pp = cp;
    }
    ndims++;
    dimchunkspecs.ndims = ndims;
    dimchunkspecs.idimids = (int *) emalloc(ndims * sizeof(int));
    dimchunkspecs.chunksizes = (size_t *) emalloc(ndims * sizeof(size_t));
    /* Look up dimension ids and assign chunksizes */
    pp = spec;
    np = spec;
    idim = 0;
    for(cp = spec; ; cp++) {
	if(*cp == '\0' || (*cp == ',' && *pp != '\\')) { /* found end of "dim/nn" part */
	    char* dimname = 0;
	    char *dp;
	    int dimid;
	    size_t chunksize;
	 
	    for(; pp > np && *pp != '/'; pp--) { /* look backwards for "/" */
		continue;
	    }
	    if(*pp != '/') {	/* no '/' found, no chunksize specified for dimension */
		ret = NC_EINVAL;
		goto done;
	    }
	    /* extract dimension name */
	    dimname = (char *) emalloc(pp - np + 1);
	    dp = dimname;
	    while(np < pp) {
		*dp++ = *np++;
	    }
	    *dp = '\0';
	    /* look up dimension id from dimension pathname */
	    ret = nc_inq_dimid2(igrp, dimname, &dimid);
	    if(ret != NC_NOERR)
		{if(dimname) free(dimname); goto done;}
	    dimchunkspecs.idimids[idim] = dimid;
	    /* parse and assign corresponding chunksize */
	    pp++; /* now points to first digit of chunksize, ',', or '\0' */
	    if(*pp == ',' || *pp == '\0') { /* no size specified, use dim len */
		size_t dimlen;
		ret = nc_inq_dimlen(igrp, dimid, &dimlen);
		if(ret != NC_NOERR)
		    {if(dimname) free(dimname); goto done;}
		chunksize = dimlen;
	    } else {	      /* convert nnn string to long long integer */
		char *ep;
#ifdef HAVE_STRTOLL
		long long val = strtoll(pp, &ep, 0);
#else
		long long val = strtol(pp, &ep, 0);
#endif
		if(ep == pp || errno == ERANGE || val < 1) /* allow chunksize bigger than dimlen */
		    {if(dimname) free(dimname); ret = NC_EINVAL; goto done;}
		chunksize = (size_t)val;
	    }
	    dimchunkspecs.chunksizes[idim] = chunksize;
	    idim++;
	    if(dimname) free(dimname);
	    dimname = NULL;
	    if(*cp == '\0')
		break;
	    /* set np to point to first char after comma */
	    np = cp + 1;
	}
	pp = cp;
    };
done:
    return ret;
}

/* Return size in chunkspec string specified for dimension corresponding to dimid, 0 if not found */
size_t
dimchunkspec_size(int indimid) {
    int idim;
    for(idim = 0; idim < dimchunkspecs.ndims; idim++) {
	if(indimid == dimchunkspecs.idimids[idim]) {
	    return dimchunkspecs.chunksizes[idim];
	}	
    }
    return 0;
}

/* Return number of dimensions for which chunking was specified in
 * chunkspec string on command line, 0 if no chunkspec string was
 * specified. */
int
dimchunkspec_ndims(void) {
    return dimchunkspecs.ndims;
}

/* Return whether chunking should be omitted, due to explicit
 * command-line specification. */
bool_t
dimchunkspec_omit(void) {
    return dimchunkspecs.omit;
}


/*
 * Parse per-variable chunkspec string and convert into varchunkspec structure.
 *   ncid: location ID of open netCDF file or group in an open file
 *   spec: string of form 
 *           var:n1,n2,...nk
 *
 *         specifying chunk size (ni) to be used for ith dimension of
 *         variable named var. Variable names may be absolute.
 *         e.g. "/grp_a/grp_a1/var".
 *         If no chunk sizes are specified, then the variable is not chunked at all.
 *
 * Returns NC_NOERR if no error, NC_EINVAL if spec has consecutive
 * unescaped commas or no chunksize specified for dimension.
 */
static int
varchunkspec_parse(int igrp, const char *spec0)
{
    int ret = NC_NOERR;
    int rank;
    int i;
    int dimids[NC_MAX_VAR_DIMS];
    struct VarChunkSpec* chunkspec = NULL;
    char* spec = NULL;
    char* p, *q; /* for walking strings */

    /* Copy spec so we can modify in place */
    spec = strdup(spec0);
    if(spec == NULL) {ret = NC_ENOMEM; goto done;}

    chunkspec = calloc(1,sizeof(struct VarChunkSpec));
    if(chunkspec == NULL) {ret = NC_ENOMEM; goto done;}

    chunkspec->igrpid = igrp;

    /* First, find the end of the variable part */
    p = strchr(spec,':');
    if(p == NULL)
	{ret = NC_EINVAL; goto done;}
    *p++ = '\0';

    /* Lookup the variable by name */
    ret = nc_inq_varid2(igrp, spec, &chunkspec->ivarid, &chunkspec->igrpid);
    if(ret != NC_NOERR) goto done;
    
    if(*p == '\0') {/* we have -c var: => do not chunk var */
	chunkspec->omit = 1;
        /* add the chunkspec to our list */
        listpush(varchunkspecs,chunkspec);
        chunkspec = NULL;
	goto done;
    }

    /* Iterate over dimension sizes */
    while(*p) {
	unsigned long dimsize;
	q = strchr(p,',');
	if(q == NULL) 
	    q = p + strlen(p); /* Fake the endpoint */
	else
	    *q++ = '\0';
	
	/* Scan as unsigned long */
	if(sscanf(p,"%lu",&dimsize) != 1)
	    {ret = NC_EINVAL; goto done;} /* Apparently not a valid dimension size */
	if(chunkspec->rank >= NC_MAX_VAR_DIMS) {ret = NC_EINVAL; goto done;} /* to many chunks */
	chunkspec->chunksizes[chunkspec->rank] = (size_t)dimsize;		
	chunkspec->rank++;
	p = q;
    }
    /* Now do some validity checking */
    /* Get some info about the var (from input) */
    ret = nc_inq_var(chunkspec->igrpid,chunkspec->ivarid,NULL,NULL,&rank,dimids,NULL);
    if(ret != NC_NOERR) goto done;
    
    /* 1. check # chunksizes == rank of variable */
    if(rank != chunkspec->rank) {ret = NC_EINVAL; goto done;}

    /* 2. check that chunksizes are legal for the given dimension sizes */
    for(i=0;i<rank;i++) {
	size_t len;
	ret = nc_inq_dimlen(igrp,dimids[i],&len);
	if(ret != NC_NOERR) goto done;
	if(chunkspec->chunksizes[i] > len) {ret = NC_EBADCHUNK; goto done;}
    }
    
    /* add the chunkspec to our list */
    listpush(varchunkspecs,chunkspec);
    chunkspec = NULL;

done:
    if(chunkspec != NULL)
	free(chunkspec);
    if(spec != NULL)
	free(spec);
    return ret;
}

/* Accessors */

bool_t
varchunkspec_exists(int igrpid, int ivarid)
{
    int i;
    for(i=0;i<listlength(varchunkspecs);i++) {
	struct VarChunkSpec* spec = listget(varchunkspecs,i);
	if(spec->igrpid == igrpid && spec->ivarid == ivarid)
	    return true;
    }
    return false;
}

bool_t
varchunkspec_omit(int igrpid, int ivarid)
{
    int i;
    for(i=0;i<listlength(varchunkspecs);i++) {
	struct VarChunkSpec* spec = listget(varchunkspecs,i);
	if(spec->igrpid == igrpid && spec->ivarid == ivarid)
	    return spec->omit;
    }
    return dimchunkspecs.omit;
}

size_t*
varchunkspec_chunksizes(int igrpid, int ivarid)
{
    int i;
    for(i=0;i<listlength(varchunkspecs);i++) {
	struct VarChunkSpec* spec = listget(varchunkspecs,i);
	if(spec->igrpid == igrpid && spec->ivarid == ivarid)
	    return spec->chunksizes;
    }
    return NULL;
}

size_t
varchunkspec_rank(int igrpid, int ivarid)
{
    int i;
    for(i=0;i<listlength(varchunkspecs);i++) {
	struct VarChunkSpec* spec = listget(varchunkspecs,i);
	if(spec->igrpid == igrpid && spec->ivarid == ivarid)
	    return spec->rank;
    }
    return 0;
}


