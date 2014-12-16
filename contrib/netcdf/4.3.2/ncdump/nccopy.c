/*********************************************************************
 *   Copyright 2010, University Corporation for Atmospheric Research
 *   See netcdf/README file for copying and redistribution conditions.
 *   Thanks to Philippe Poilbarbe and Antonio S. Cofiño for 
 *   compression additions.
 *   $Id: nccopy.c 400 2010-08-27 21:02:52Z russ $
 *********************************************************************/

#include "config.h"		/* for USE_NETCDF4 macro */
#include <stdlib.h>
#ifdef HAVE_GETOPT_H
#include <getopt.h>
#endif
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#include <string.h>
#include <netcdf.h>
#include "nciter.h"
#include "utils.h"
#include "chunkspec.h"
#include "dimmap.h"
#include "nccomps.h"

#ifdef _MSC_VER
#include "XGetopt.h"
#define snprintf _snprintf
int opterr;
int optind;
#endif

/* default bytes of memory we are willing to allocate for variable
 * values during copy */
#define COPY_BUFFER_SIZE (5000000)
#define COPY_CHUNKCACHE_PREEMPTION (1.0f) /* for copying, can eject fully read chunks */
#define SAME_AS_INPUT (-1)	/* default, if kind not specified */
#define CHUNK_THRESHOLD (8192)	/* non-record variables with fewer bytes don't get chunked */

#ifndef USE_NETCDF4
#define NC_CLASSIC_MODEL 0x0100 /* Enforce classic model if netCDF-4 not available. */
#endif

/* Global variables for command-line requests */
char *progname;	       /* for error messages */
static int option_kind = SAME_AS_INPUT;
static int option_deflate_level = -1;	/* default, compress output only if input compressed */
static int option_shuffle_vars = NC_NOSHUFFLE; /* default, no shuffling on compression */
static int option_fix_unlimdims = 0; /* default, preserve unlimited dimensions */
static char* option_chunkspec = 0;   /* default, no chunk specification */
static size_t option_copy_buffer_size = COPY_BUFFER_SIZE;
static size_t option_chunk_cache_size = CHUNK_CACHE_SIZE; /* default from config.h */
static size_t option_chunk_cache_nelems = CHUNK_CACHE_NELEMS; /* default from config.h */
static int option_read_diskless = 0; /* default, don't read input into memory on open */
static int option_write_diskless = 0; /* default, don't write output to diskless file */
static int option_min_chunk_bytes = CHUNK_THRESHOLD; /* default, don't chunk variable if prod of
						      * chunksizes of its dimensions is smaller
						      * than this */
static int option_nlgrps = 0;		    /* Number of groups specified with -g
					     * option on command line */
static char** option_lgrps = 0;		    /* list of group names specified with -g
					     * option on command line */
static idnode_t* option_grpids = 0; /* list of grpids matching list specified with -g option */
static bool_t option_grpstruct = false; /* if -g set, copy structure for non-selected groups */
static int option_nlvars = 0; /* Number of variables specified with -v * option on command line */
static char** option_lvars = 0;         /* list of variable names specified with -v
                                         * option on command line */
static bool_t option_varstruct = false;   /* if -v set, copy structure for non-selected vars */
static int option_compute_chunkcaches = 0; /* default, don't try still flaky estimate of
					    * chunk cache for each variable */

/* get group id in output corresponding to group igrp in input,
 * given parent group id (or root group id) parid in output. */
static int
get_grpid(int igrp, int parid, int *ogrpp) {
    int stat = NC_NOERR;
    int ogid = parid;		/* like igrp but in output file */
#ifdef USE_NETCDF4
    int inparid;

    /* if not root group, get corresponding output groupid from group name */
    stat = nc_inq_grp_parent(igrp, &inparid);
    if(stat == NC_NOERR) {	/* not root group */
	char grpname[NC_MAX_NAME + 1];
	NC_CHECK(nc_inq_grpname(igrp, grpname));
	NC_CHECK(nc_inq_grp_ncid(parid, grpname, &ogid));
    } else if(stat == NC_ENOGRP) { /* root group */
	stat = NC_NOERR;
    } else {
	NC_CHECK(stat);
    }
#endif	/* USE_NETCDF4 */
    *ogrpp = ogid;
    return stat;
}

/* Return size in bytes of a variable value */
static size_t
val_size(int grpid, int varid) {
    nc_type vartype;
    size_t value_size;
    NC_CHECK(nc_inq_vartype(grpid, varid, &vartype));
    NC_CHECK(nc_inq_type(grpid, vartype, NULL, &value_size));
    return value_size;
}

#ifdef USE_NETCDF4
/* Get parent id needed to define a new group from its full name in an
 * open file identified by ncid.  Assumes all intermediate groups are
 * already defined.  */
static int
nc_inq_parid(int ncid, const char *fullname, int *locidp) {
    char *parent = strdup(fullname);
    char *slash = "/";		/* groupname separator */
    char *last_slash;
    if(parent == NULL) {
	return NC_ENOMEM;	/* exits */
    }
    last_slash = strrchr(parent, '/');
    if(last_slash == parent || last_slash == NULL) {	/* parent is root */
	free(parent);
	parent = strdup(slash);
    } else {
	*last_slash = '\0';	/* truncate to get parent name */
    }
    NC_CHECK(nc_inq_grp_full_ncid(ncid, parent, locidp));
       free(parent);
    return NC_NOERR;
}

/* Return size of chunk in bytes for a variable varid in a group igrp, or 0 if
 * layout is contiguous */
static int
inq_var_chunksize(int igrp, int varid, size_t* chunksizep) {
    int stat = NC_NOERR;
    int ndims;
    size_t *chunksizes;
    int dim;
    int contig = 1;
    nc_type vartype;
    size_t value_size;
    size_t prod;

    NC_CHECK(nc_inq_vartype(igrp, varid, &vartype));
    /* from type, get size in memory needed for each value */
    NC_CHECK(nc_inq_type(igrp, vartype, NULL, &value_size));
    prod = value_size;
    NC_CHECK(nc_inq_varndims(igrp, varid, &ndims));
    chunksizes = (size_t *) emalloc((ndims + 1) * sizeof(size_t));
    if(ndims > 0) {
	NC_CHECK(nc_inq_var_chunking(igrp, varid, &contig, NULL));
    }
    if(contig == 1) {
	*chunksizep = 0;
    } else {
	NC_CHECK(nc_inq_var_chunking(igrp, varid, &contig, chunksizes));
	for(dim = 0; dim < ndims; dim++) {
	    prod *= chunksizes[dim];
	}
	*chunksizep = prod;
    }
    free(chunksizes);
    return stat;
}

/* Return estimated number of elems required in chunk cache and
 * estimated size of chunk cache adequate to efficiently copy input
 * variable ivarid to output variable ovarid, which may have different
 * chunk size and shape */
static int
inq_var_chunking_params(int igrp, int ivarid, int ogrp, int ovarid,
			size_t* chunkcache_sizep,
			size_t *chunkcache_nelemsp,
                        float * chunkcache_preemptionp)
{
    int stat = NC_NOERR;
    int ndims;
    size_t *ichunksizes, *ochunksizes;
    int dim;
    int icontig = 1, ocontig = 1;
    nc_type vartype;
    size_t value_size;
    size_t prod, iprod, oprod;
    size_t nelems;
    *chunkcache_nelemsp = CHUNK_CACHE_NELEMS;
    *chunkcache_sizep = CHUNK_CACHE_SIZE;
    *chunkcache_preemptionp = COPY_CHUNKCACHE_PREEMPTION;

    NC_CHECK(nc_inq_varndims(igrp, ivarid, &ndims));
    if(ndims > 0) {
	NC_CHECK(nc_inq_var_chunking(igrp, ivarid, &icontig, NULL));
	NC_CHECK(nc_inq_var_chunking(ogrp, ovarid, &ocontig, NULL));
    }
    if(icontig == 1 && ocontig == 1) { /* no chunking in input or output */
	*chunkcache_nelemsp = 0;
	*chunkcache_sizep = 0;
	*chunkcache_preemptionp = 0;
	return stat;
    }

    NC_CHECK(nc_inq_vartype(igrp, ivarid, &vartype));
    NC_CHECK(nc_inq_type(igrp, vartype, NULL, &value_size));
    iprod = value_size;

    if(icontig == 0 && ocontig == 1) { /* chunking only in input */
	*chunkcache_nelemsp = 1;       /* read one input chunk at a time */
	*chunkcache_sizep = iprod;
	*chunkcache_preemptionp = 1.0f;
	return stat;
    }

    ichunksizes = (size_t *) emalloc((ndims + 1) * sizeof(size_t));
    if(icontig == 1) { /* if input contiguous, treat as if chunked on
			* first dimension */
	ichunksizes[0] = 1;
	for(dim = 1; dim < ndims; dim++) {
	    ichunksizes[dim] = dim;
	}
    } else {
	NC_CHECK(nc_inq_var_chunking(igrp, ivarid, &icontig, ichunksizes));
    }

    /* now can assume chunking in both input and output */
    ochunksizes = (size_t *) emalloc((ndims + 1) * sizeof(size_t));
    NC_CHECK(nc_inq_var_chunking(ogrp, ovarid, &ocontig, ochunksizes));

    nelems = 1;
    oprod = value_size;
    for(dim = 0; dim < ndims; dim++) {
	nelems += 1 + (ichunksizes[dim] - 1) / ochunksizes[dim];
	iprod *= ichunksizes[dim];
	oprod *= ochunksizes[dim];
    }
    prod = iprod + oprod * (nelems - 1);
    *chunkcache_nelemsp = nelems;
    *chunkcache_sizep = prod;
    free(ichunksizes);
    free(ochunksizes);
    return stat;
}

/* Forward declaration, because copy_type, copy_vlen_type call each other */
static int copy_type(int igrp, nc_type typeid, int ogrp);

/* 
 * copy a user-defined variable length type in the group igrp to the
 * group ogrp
 */
static int
copy_vlen_type(int igrp, nc_type itype, int ogrp)
{
    int stat = NC_NOERR; 
    nc_type ibasetype;
    nc_type obasetype;		/* base type in target group */
    char name[NC_MAX_NAME];
    size_t size;
    char basename[NC_MAX_NAME];
    size_t basesize;
    nc_type vlen_type;

    NC_CHECK(nc_inq_vlen(igrp, itype, name, &size, &ibasetype));
    /* to get base type id in target group, use name of base type in
     * source group */
    NC_CHECK(nc_inq_type(igrp, ibasetype, basename, &basesize));
    stat = nc_inq_typeid(ogrp, basename, &obasetype);
    /* if no such type, create it now */
    if(stat == NC_EBADTYPE) {
	NC_CHECK(copy_type(igrp, ibasetype, ogrp));
	stat = nc_inq_typeid(ogrp, basename, &obasetype);
    }
    NC_CHECK(stat);

    /* Now we know base type exists in output and we know its type id */
    NC_CHECK(nc_def_vlen(ogrp, name, obasetype, &vlen_type));

    return stat;
}

/* 
 * copy a user-defined opaque type in the group igrp to the group ogrp
 */
static int
copy_opaque_type(int igrp, nc_type itype, int ogrp)
{
    int stat = NC_NOERR; 
    nc_type otype;
    char name[NC_MAX_NAME];
    size_t size;

    NC_CHECK(nc_inq_opaque(igrp, itype, name, &size));
    NC_CHECK(nc_def_opaque(ogrp, size, name, &otype));

    return stat;
}

/* 
 * copy a user-defined enum type in the group igrp to the group ogrp
 */
static int
copy_enum_type(int igrp, nc_type itype, int ogrp)
{
    int stat = NC_NOERR; 
    nc_type otype;
    nc_type basetype;
    size_t basesize;
    size_t nmembers;
    char name[NC_MAX_NAME];
    int i;

    NC_CHECK(nc_inq_enum(igrp, itype, name, &basetype, &basesize, &nmembers));
    NC_CHECK(nc_def_enum(ogrp, basetype, name, &otype));
    for(i = 0; i < nmembers; i++) { /* insert enum members */
	char ename[NC_MAX_NAME];
	long long val;		/* large enough to hold any integer type */
	NC_CHECK(nc_inq_enum_member(igrp, itype, i, ename, &val));
	NC_CHECK(nc_insert_enum(ogrp, otype, ename, &val));
    }
    return stat;
}

/* 
 * copy a user-defined compound type in the group igrp to the group ogrp
 */
static int
copy_compound_type(int igrp, nc_type itype, int ogrp)
{
    int stat = NC_NOERR; 
    char name[NC_MAX_NAME];
    size_t size;
    size_t nfields;
    nc_type otype;
    int fid;

    NC_CHECK(nc_inq_compound(igrp, itype, name, &size, &nfields));
    NC_CHECK(nc_def_compound(ogrp, size, name, &otype));

    for (fid = 0; fid < nfields; fid++) {
	char fname[NC_MAX_NAME];
	char ftypename[NC_MAX_NAME];
	size_t foff;
	nc_type iftype, oftype;
	int fndims;

	NC_CHECK(nc_inq_compound_field(igrp, itype, fid, fname, &foff, &iftype, &fndims, NULL));
	/* type ids in source don't necessarily correspond to same
	 * typeids in destination, so look up destination typeid by using
	 * field type name */
	NC_CHECK(nc_inq_type(igrp, iftype, ftypename, NULL));
	NC_CHECK(nc_inq_typeid(ogrp, ftypename, &oftype));
	if(fndims == 0) {
	    NC_CHECK(nc_insert_compound(ogrp, otype, fname, foff, oftype));
	} else {		/* field is array type */
	    int *fdimsizes;
	    fdimsizes = (int *) emalloc((fndims + 1) * sizeof(int));
	    stat = nc_inq_compound_field(igrp, itype, fid, NULL, NULL, NULL, 
					 NULL, fdimsizes);
	    NC_CHECK(nc_insert_array_compound(ogrp, otype, fname, foff, oftype, fndims, fdimsizes));
	    free(fdimsizes);
	}
    }
    return stat;
}


/* 
 * copy a user-defined type in the group igrp to the group ogrp
 */
static int
copy_type(int igrp, nc_type typeid, int ogrp)
{
    int stat = NC_NOERR; 
    nc_type type_class;

    NC_CHECK(nc_inq_user_type(igrp, typeid, NULL, NULL, NULL, NULL, &type_class));

    switch(type_class) {
    case NC_VLEN:
	NC_CHECK(copy_vlen_type(igrp, typeid, ogrp));
	break;
    case NC_OPAQUE:
	NC_CHECK(copy_opaque_type(igrp, typeid, ogrp));
	break;
    case NC_ENUM:
	NC_CHECK(copy_enum_type(igrp, typeid, ogrp));
	break;
    case NC_COMPOUND:
	NC_CHECK(copy_compound_type(igrp, typeid, ogrp));
	break;
    default:
	NC_CHECK(NC_EBADTYPE);
    }
    return stat;
}

/* Copy a group and all its subgroups, recursively, from iroot to
 * oroot, the ncids of input file and output file.  This just creates
 * all the groups in the destination, but doesn't copy anything that's
 * in the groups yet. */
static int
copy_groups(int iroot, int oroot)
{
    int stat = NC_NOERR;
    int numgrps;
    int *grpids;
    int i;

    /* get total number of groups and their ids, including all descendants */
    NC_CHECK(nc_inq_grps_full(iroot, &numgrps, NULL));
    if(numgrps > 1) {		/* there's always 1 root group */
	grpids = emalloc(numgrps * sizeof(int));
	NC_CHECK(nc_inq_grps_full(iroot, NULL, grpids));
	/* create corresponding new groups in ogrp, except for root group */
	for(i = 1; i < numgrps; i++) {
	    char *grpname_full;
	    char grpname[NC_MAX_NAME];
	    size_t len_name;
	    int ogid = 0, oparid = 0, iparid = 0;
	    /* get full group name of input group */
	    NC_CHECK(nc_inq_grpname(grpids[i], grpname));
	    if (option_grpstruct || group_wanted(grpids[i], option_nlgrps, option_grpids)) {
	        NC_CHECK(nc_inq_grpname_full(grpids[i], &len_name, NULL));
		grpname_full = emalloc(len_name + 1);
		NC_CHECK(nc_inq_grpname_full(grpids[i], &len_name, grpname_full));
		/* Make sure, the parent group is also wanted (root group is always wanted) */
		NC_CHECK(nc_inq_parid(iroot, grpname_full, &iparid));
		if (!option_grpstruct && !group_wanted(iparid, option_nlgrps, option_grpids) 
		    && iparid != iroot) {
		    error("ERROR: trying to copy a group but not the parent: %s", grpname_full);
		}
		/* get id of parent group of corresponding group in output.
		 * Note that this exists, because nc_inq_groups returned
		 * grpids in preorder, so parents are always copied before
		 * their subgroups */
		NC_CHECK(nc_inq_parid(oroot, grpname_full, &oparid));
		NC_CHECK(nc_inq_grpname(grpids[i], grpname));
		/* define corresponding group in output */
		NC_CHECK(nc_def_grp(oparid, grpname, &ogid));
		free(grpname_full);
	    }
	}
	free(grpids);
    }
    return stat;    
}

/* 
 * Copy the user-defined types in this group (igrp) and all its
 * subgroups, recursively, to corresponding group in output (ogrp)
 */
static int
copy_types(int igrp, int ogrp)
{
    int stat = NC_NOERR; 
    int ntypes;
    nc_type *types = NULL;
    int numgrps;
    int *grpids = NULL;
    int i;

    NC_CHECK(nc_inq_typeids(igrp, &ntypes, NULL));

    if(ntypes > 0) {
	types = (nc_type *) emalloc(ntypes * sizeof(nc_type));
	NC_CHECK(nc_inq_typeids(igrp, &ntypes, types));
	for (i = 0; i < ntypes; i++) {
	    NC_CHECK(copy_type(igrp, types[i], ogrp));
	}
	free(types);
    }

    /* Copy types from subgroups */
    NC_CHECK(nc_inq_grps(igrp, &numgrps, NULL));
    if(numgrps > 0) {
	grpids = (int *)emalloc(sizeof(int) * numgrps);
	NC_CHECK(nc_inq_grps(igrp, &numgrps, grpids));
	for(i = 0; i < numgrps; i++) {
	    if (option_grpstruct || group_wanted(grpids[i], option_nlgrps, option_grpids)) {
		int ogid;
		/* get groupid in output corresponding to grpids[i] in
		 * input, given parent group (or root group) ogrp in
		 * output */
		NC_CHECK(get_grpid(grpids[i], ogrp, &ogid));
		NC_CHECK(copy_types(grpids[i], ogid));
	    }
	}
	free(grpids);
    }
    return stat;
}

/* Copy all netCDF-4 specific variable properties such as chunking,
 * endianness, deflation, checksumming, fill, etc. */
static int
copy_var_specials(int igrp, int varid, int ogrp, int o_varid)
{
    int stat = NC_NOERR;
    {				/* handle chunking parameters */
	int ndims;
	NC_CHECK(nc_inq_varndims(igrp, varid, &ndims));
	if (ndims > 0) {		/* no chunking for scalar variables */
	    int contig = 0;
	    size_t *chunkp = (size_t *) emalloc(ndims * sizeof(size_t));
	    int *dimids = (int *) emalloc(ndims * sizeof(int));
	    int idim;
	     /* size of a chunk: product of dimension chunksizes and size of value */ 
	    size_t csprod = val_size(ogrp, o_varid);
	    int is_unlimited = 0;
	    NC_CHECK(nc_inq_var_chunking(igrp, varid, &contig, chunkp));
	    NC_CHECK(nc_inq_vardimid(igrp, varid, dimids));

	    for(idim = 0; idim < ndims; idim++) {
		int idimid = dimids[idim];
		int odimid = dimmap_odimid(idimid);
		size_t chunksize = chunkspec_size(idimid);
		if(chunksize > 0) { /* found in chunkspec */
		    chunkp[idim] = chunksize;
		}
		csprod *= chunkp[idim];
		if(dimmap_ounlim(odimid))
		    is_unlimited = 1;
	    }
	    /* Explicitly set chunking, even if default */
	    /* If product of chunksizes is too small and no unlimited
	     * dimensions used, don't chunk.  Also if chunking
	     * explicitly turned off with chunk spec, don't chunk. */
	    if ((csprod < option_min_chunk_bytes && !is_unlimited) || contig == 1
		|| chunkspec_omit() == true) {
		NC_CHECK(nc_def_var_chunking(ogrp, o_varid, NC_CONTIGUOUS, NULL));
	    } else {
		NC_CHECK(nc_def_var_chunking(ogrp, o_varid, NC_CHUNKED, chunkp));
	    }
	    free(dimids);
	    free(chunkp);
	}
    }
    { /* handle compression parameters, copying from input, overriding
       * with command-line options */
	int shuffle_in=0, deflate_in=0, deflate_level_in=0;
	int shuffle_out=0, deflate_out=0, deflate_level_out=0;
	if(option_deflate_level != 0) {
	    NC_CHECK(nc_inq_var_deflate(igrp, varid, &shuffle_in, &deflate_in, &deflate_level_in));
	    if(option_deflate_level == -1) { /* not specified, copy input compression and shuffling */
		shuffle_out = shuffle_in;
		deflate_out = deflate_in;
		deflate_level_out = deflate_level_in;
	    } else if(option_deflate_level > 0) { /* change to specified compression, shuffling */
		shuffle_out = option_shuffle_vars;
		deflate_out=1;
		deflate_level_out = option_deflate_level;
	    }
	    NC_CHECK(nc_def_var_deflate(ogrp, o_varid, shuffle_out, deflate_out, deflate_level_out));
	}
    }
    {				/* handle checksum parameters */
	int fletcher32 = 0;
	NC_CHECK(nc_inq_var_fletcher32(igrp, varid, &fletcher32));
	if(fletcher32 != 0) {
	    NC_CHECK(nc_def_var_fletcher32(ogrp, o_varid, fletcher32));
	}
    }
    {				/* handle endianness */
	int endianness = 0;
	NC_CHECK(nc_inq_var_endian(igrp, varid, &endianness));
	if(endianness != NC_ENDIAN_NATIVE) { /* native is the default */
	    NC_CHECK(nc_def_var_endian(ogrp, o_varid, endianness));
	}
    }
    return stat;
}

/* Set output variable o_varid (in group ogrp) to use chunking
 * specified on command line, only called for classic format input and
 * netCDF-4 format output, so no existing chunk lengths to override. */
static int
set_var_chunked(int ogrp, int o_varid)
{
    int stat = NC_NOERR;
    int ndims;
    int odim;
    size_t chunk_threshold = CHUNK_THRESHOLD;

    if(chunkspec_ndims() == 0) 	/* no chunking specified on command line */
	return stat;
    NC_CHECK(nc_inq_varndims(ogrp, o_varid, &ndims));

    if (ndims > 0) {		/* no chunking for scalar variables */
	int chunked = 0;
	int *dimids = (int *) emalloc(ndims * sizeof(int));
	size_t varsize;
	nc_type vartype;
	size_t value_size;
	int is_unlimited = 0;

	NC_CHECK(nc_inq_vardimid (ogrp, o_varid, dimids));
	NC_CHECK(nc_inq_vartype(ogrp, o_varid, &vartype));
	/* from type, get size in memory needed for each value */
	NC_CHECK(nc_inq_type(ogrp, vartype, NULL, &value_size));
	varsize = value_size;

	/* Determine if this variable should be chunked.  A variable
	 * should be chunked if any of its dims are in command-line
	 * chunk spec. It will also be chunked if any of its
	 * dims are unlimited. */
	for(odim = 0; odim < ndims; odim++) {
	    int odimid = dimids[odim];
	    int idimid = dimmap_idimid(odimid); /* corresponding dimid in input file */
	    if(dimmap_ounlim(odimid))
		is_unlimited = 1;
	    if(idimid != -1) {
		size_t chunksize = chunkspec_size(idimid); /* from chunkspec */
		size_t dimlen;
		NC_CHECK(nc_inq_dimlen(ogrp, odimid, &dimlen));
		if( (chunksize > 0) || dimmap_ounlim(odimid)) {
		    chunked = 1;		    
		}
		varsize *= dimlen;
	    }
	}
	/* Don't chunk small variables that don't use an unlimited
	 * dimension. */
	if(varsize < chunk_threshold && !is_unlimited)
	    chunked = 0;

	if(chunked) {
	    /* Allocate chunksizes and set defaults to dimsize for any
	     * dimensions not mentioned in chunkspec. */
	    size_t *chunkp = (size_t *) emalloc(ndims * sizeof(size_t));
	    for(odim = 0; odim < ndims; odim++) {
		int odimid = dimids[odim];
		int idimid = dimmap_idimid(odimid);
		size_t chunksize = chunkspec_size(idimid);
		if(chunksize > 0) {
		    chunkp[odim] = chunksize;
		} else {
		    NC_CHECK(nc_inq_dimlen(ogrp, odimid, &chunkp[odim]));
		}
	    }
	    NC_CHECK(nc_def_var_chunking(ogrp, o_varid, NC_CHUNKED, chunkp));
	    free(chunkp);
	}
	free(dimids);
    }
    return stat;
}

/* Set variable to compression specified on command line */
static int
set_var_compressed(int ogrp, int o_varid)
{
    int stat = NC_NOERR;
    if (option_deflate_level > 0) {
	int deflate = 1;
	NC_CHECK(nc_def_var_deflate(ogrp, o_varid, option_shuffle_vars, deflate, option_deflate_level));
    }
    return stat;
}

/* Release the variable chunk cache allocated for variable varid in
 * group grp.  This is not necessary, but will save some memory when
 * processing one variable at a time.  */
#ifdef UNUSED
static int
free_var_chunk_cache(int grp, int varid)
{
    int stat = NC_NOERR;
    size_t chunk_cache_size = 1;
    size_t cache_nelems = 1;
    float cache_preemp = 0;
    int kind;
    NC_CHECK(nc_inq_format(grp, &kind));
    if(kind == NC_FORMAT_NETCDF4 || kind == NC_FORMAT_NETCDF4_CLASSIC) {
	int contig = 1;
	NC_CHECK(nc_inq_var_chunking(grp, varid, &contig, NULL));
	if(contig == 0) {	/* chunked */
	    NC_CHECK(nc_set_var_chunk_cache(grp, varid, chunk_cache_size, cache_nelems, cache_preemp));
	}
    }
    return stat;
}
#endif

#endif /* USE_NETCDF4 */

/* Copy dimensions from group igrp to group ogrp, also associate input
 * dimids with output dimids (they need not match, because the input
 * dimensions may have been defined in a different order than we define
 * the output dimensions here. */
static int
copy_dims(int igrp, int ogrp)
{
    int stat = NC_NOERR;
    int ndims;
    int dgrp;
#ifdef USE_NETCDF4
    int nunlims;
    int *dimids;
    int *unlimids;
#else
    int unlimid;
#endif /* USE_NETCDF4 */    

    NC_CHECK(nc_inq_ndims(igrp, &ndims));

#ifdef USE_NETCDF4
   /* In netCDF-4 files, dimids may not be sequential because they
    * may be defined in various groups, and we are only looking at one
    * group at a time. */
    /* Find the dimension ids in this group, don't include parents. */
    dimids = (int *) emalloc((ndims + 1) * sizeof(int));
    NC_CHECK(nc_inq_dimids(igrp, NULL, dimids, 0));
    /* Find the number of unlimited dimensions and get their IDs */
    NC_CHECK(nc_inq_unlimdims(igrp, &nunlims, NULL));
    unlimids = (int *) emalloc((nunlims + 1) * sizeof(int));
    NC_CHECK(nc_inq_unlimdims(igrp, NULL, unlimids));
#else
    NC_CHECK(nc_inq_unlimdim(igrp, &unlimid));
#endif /* USE_NETCDF4 */

    /* Copy each dimension to output, including unlimited dimension(s) */
    for (dgrp = 0; dgrp < ndims; dgrp++) {
	char name[NC_MAX_NAME];
	size_t length;
	int i_is_unlim;
	int o_is_unlim;
	int idimid, odimid;
#ifdef USE_NETCDF4
	int uld;
#endif

	i_is_unlim = 0;
#ifdef USE_NETCDF4
	idimid = dimids[dgrp];
	for (uld = 0; uld < nunlims; uld++) {
	    if(idimid == unlimids[uld]) {
		i_is_unlim = 1;
		break;
	    }	  
	}
#else
	idimid = dgrp;
	if(unlimid != -1 && (idimid == unlimid)) {
	    i_is_unlim = 1;
	}
#endif /* USE_NETCDF4 */

	stat = nc_inq_dim(igrp, idimid, name, &length);
	if (stat == NC_EDIMSIZE && sizeof(size_t) < 8) {
	    error("dimension \"%s\" requires 64-bit platform", name);
	}	
	NC_CHECK(stat);
	o_is_unlim = i_is_unlim;
	if(i_is_unlim && !option_fix_unlimdims) {
	    NC_CHECK(nc_def_dim(ogrp, name, NC_UNLIMITED, &odimid));
	} else {
	    NC_CHECK(nc_def_dim(ogrp, name, length, &odimid));
	    o_is_unlim = 0;
	}
	/* Store (idimid, odimid) mapping for later use, also whether unlimited */
	dimmap_store(idimid, odimid, i_is_unlim, o_is_unlim);
    }
#ifdef USE_NETCDF4
    free(dimids);
    free(unlimids);
#endif /* USE_NETCDF4 */    
    return stat;
}

/* Copy the attributes for variable ivar in group igrp to variable
 * ovar in group ogrp.  Global (group) attributes are specified by
 * using the varid NC_GLOBAL */
static int
copy_atts(int igrp, int ivar, int ogrp, int ovar)
{
    int natts;
    int iatt;
    int stat = NC_NOERR;

    NC_CHECK(nc_inq_varnatts(igrp, ivar, &natts));
    
    for(iatt = 0; iatt < natts; iatt++) {
	char name[NC_MAX_NAME];
	NC_CHECK(nc_inq_attname(igrp, ivar, iatt, name));
	NC_CHECK(nc_copy_att(igrp, ivar, name, ogrp, ovar));
    }
    return stat;
}

/* copy the schema for a single variable in group igrp to group ogrp */
static int
copy_var(int igrp, int varid, int ogrp)
{
    int stat = NC_NOERR;
    int ndims;
    int *idimids;		/* ids of dims for input variable */
    int *odimids;		/* ids of dims for output variable */
    char name[NC_MAX_NAME];
    nc_type typeid, o_typeid;
    int natts;
    int i;
    int o_varid;

    NC_CHECK(nc_inq_varndims(igrp, varid, &ndims));
    idimids = (int *) emalloc((ndims + 1) * sizeof(int));
    NC_CHECK(nc_inq_var(igrp, varid, name, &typeid, NULL, idimids, &natts));
    o_typeid = typeid;
#ifdef USE_NETCDF4
    if (typeid > NC_STRING) {	/* user-defined type */
	/* type ids in source don't necessarily correspond to same
	 * typeids in destination, so look up destination typeid by
	 * using type name */
	char type_name[NC_MAX_NAME];
	NC_CHECK(nc_inq_type(igrp, typeid, type_name, NULL));
	NC_CHECK(nc_inq_typeid(ogrp, type_name, &o_typeid));
    }
#endif	/* USE_NETCDF4 */

    /* get the corresponding dimids in the output file */
    odimids = (int *) emalloc((ndims + 1) * sizeof(int));
    for(i = 0; i < ndims; i++) {
	odimids[i] = dimmap_odimid(idimids[i]);
	if(odimids[i] == -1) {
	    error("Oops, no dimension in output associated with input dimid %d", idimids[i]);
	}
    }

    /* define the output variable */
    NC_CHECK(nc_def_var(ogrp, name, o_typeid, ndims, odimids, &o_varid));
    /* attach the variable attributes to the output variable */
    NC_CHECK(copy_atts(igrp, varid, ogrp, o_varid));
#ifdef USE_NETCDF4    
    {
	int inkind;
	int outkind;
	NC_CHECK(nc_inq_format(igrp, &inkind));
	NC_CHECK(nc_inq_format(ogrp, &outkind));
	if(outkind == NC_FORMAT_NETCDF4 || outkind == NC_FORMAT_NETCDF4_CLASSIC) {
	    if((inkind == NC_FORMAT_NETCDF4 || inkind == NC_FORMAT_NETCDF4_CLASSIC)) {
		/* Copy all netCDF-4 specific variable properties such as
		 * chunking, endianness, deflation, checksumming, fill, etc. */
		NC_CHECK(copy_var_specials(igrp, varid, ogrp, o_varid));
	    } else {
		/* Set chunking if specified in command line option */
		NC_CHECK(set_var_chunked(ogrp, o_varid));
		/* Set compression if specified in command line option */
		NC_CHECK(set_var_compressed(ogrp, o_varid));
	    }
	}
    }
#endif	/* USE_NETCDF4 */
    free(idimids);
    free(odimids);
    return stat;
}

/* copy the schema for all the variables in group igrp to group ogrp */
static int
copy_vars(int igrp, int ogrp)
{
    int stat = NC_NOERR;
    int nvars;
    int varid;

    int iv;			/* variable number */
    idnode_t* vlist = 0;		/* list for vars specified with -v option */

    /*
     * If any vars were specified with -v option, get list of
     * associated variable ids relative to this group.  Assume vars
     * specified with syntax like "grp1/grp2/varname" or
     * "/grp1/grp2/varname" if they are in groups.
     */
    vlist = newidlist();	/* list for vars specified with -v option */
    for (iv=0; iv < option_nlvars; iv++) {
        if(nc_inq_gvarid(igrp, option_lvars[iv], &varid) == NC_NOERR)
            idadd(vlist, varid);
    }
    
    NC_CHECK(nc_inq_nvars(igrp, &nvars));
    for (varid = 0; varid < nvars; varid++) {
	if (!option_varstruct && option_nlvars > 0 && ! idmember(vlist, varid))
            continue;
	NC_CHECK(copy_var(igrp, varid, ogrp));
    }
    freeidlist(vlist);
    return stat;
}

/* Copy the schema in a group and all its subgroups, recursively, from
 * group igrp in input to parent group ogrp in destination.  Use
 * dimmap array to map input dimids to output dimids. */
static int
copy_schema(int igrp, int ogrp) 
{
    int stat = NC_NOERR;
    int ogid;			/* like igrp but in output file */

    /* get groupid in output corresponding to group igrp in input,
     * given parent group (or root group) ogrp in output */
    NC_CHECK(get_grpid(igrp, ogrp, &ogid));

    NC_CHECK(copy_dims(igrp, ogid));
    NC_CHECK(copy_atts(igrp, NC_GLOBAL, ogid, NC_GLOBAL));
    NC_CHECK(copy_vars(igrp, ogid));
#ifdef USE_NETCDF4    
    {
	int numgrps;
	int *grpids;
	int i;
	/* Copy schema from subgroups */
	stat = nc_inq_grps(igrp, &numgrps, NULL);
	grpids = (int *)emalloc((numgrps + 1) * sizeof(int));
	NC_CHECK(nc_inq_grps(igrp, &numgrps, grpids));
	
	for(i = 0; i < numgrps; i++) {
	    if (option_grpstruct || group_wanted(grpids[i], option_nlgrps, option_grpids)) {
	        NC_CHECK(copy_schema(grpids[i], ogid));
	    }
	}
	free(grpids);
    }
#endif	/* USE_NETCDF4 */
    return stat;    
}

/* Return number of values for a variable varid in a group igrp */
static int
inq_nvals(int igrp, int varid, long long *nvalsp) {
    int stat = NC_NOERR;
    int ndims;
    int *dimids;
    int dim;
    long long nvals = 1;

    NC_CHECK(nc_inq_varndims(igrp, varid, &ndims));
    dimids = (int *) emalloc((ndims + 1) * sizeof(int));
    NC_CHECK(nc_inq_vardimid (igrp, varid, dimids));
    for(dim = 0; dim < ndims; dim++) {
	size_t len;
	NC_CHECK(nc_inq_dimlen(igrp, dimids[dim], &len));
	nvals *= len;
    }
    if(nvalsp)
	*nvalsp = nvals;
    free(dimids);
    return stat;
}

/* Copy data from variable varid in group igrp to corresponding group
 * ogrp. */
static int
copy_var_data(int igrp, int varid, int ogrp) {
    int stat = NC_NOERR;
    nc_type vartype;
    long long nvalues;		/* number of values for this variable */
    size_t ntoget;		/* number of values to access this iteration */
    size_t value_size;		/* size of a single value of this variable */
    static void *buf = 0;	/* buffer for the variable values */
    char varname[NC_MAX_NAME];
    int ovarid;
    size_t *start;
    size_t *count;
    nciter_t *iterp;		/* opaque structure for iteration status */
    int do_realloc = 0;
#ifdef USE_NETCDF4    
    int okind;
    size_t chunksize;
#endif

    NC_CHECK(inq_nvals(igrp, varid, &nvalues));
    if(nvalues == 0)
	return stat;
    /* get corresponding output variable */
    NC_CHECK(nc_inq_varname(igrp, varid, varname));
    NC_CHECK(nc_inq_varid(ogrp, varname, &ovarid));
    NC_CHECK(nc_inq_vartype(igrp, varid, &vartype));
    value_size = val_size(igrp, varid);
    if(value_size > option_copy_buffer_size) {
	option_copy_buffer_size = value_size;
	do_realloc = 1;
    }
#ifdef USE_NETCDF4    
    NC_CHECK(nc_inq_format(ogrp, &okind));
    if(okind == NC_FORMAT_NETCDF4 || okind == NC_FORMAT_NETCDF4_CLASSIC) {
	/* if this variable chunked, set variable chunk cache size */ 
	int contig = 1;
	NC_CHECK(nc_inq_var_chunking(ogrp, ovarid, &contig, NULL));
	if(contig == 0) {	/* chunked */
	    if(option_compute_chunkcaches) {
		/* Try to estimate variable-specific chunk cache,
		 * depending on specific size and shape of this
		 * variable's chunks.  This doesn't work yet. */
		size_t chunkcache_size, chunkcache_nelems;
		float chunkcache_preemption;
		NC_CHECK(inq_var_chunking_params(igrp, varid, ogrp, ovarid,
						 &chunkcache_size, 
						 &chunkcache_nelems, 
						 &chunkcache_preemption));
		NC_CHECK(nc_set_var_chunk_cache(ogrp, ovarid, 
						chunkcache_size, 
						chunkcache_nelems, 
						chunkcache_preemption)); 
	    } else {		
		/* by default, use same chunk cache for all chunked variables */
		NC_CHECK(nc_set_var_chunk_cache(ogrp, ovarid, 
						option_chunk_cache_size,
						option_chunk_cache_nelems,
						COPY_CHUNKCACHE_PREEMPTION));
	    }
	}
    }
    /* For chunked variables, option_copy_buffer_size must also be at least as large as
     * size of a chunk in input, otherwise resize it. */
    {
	NC_CHECK(inq_var_chunksize(igrp, varid, &chunksize));
	if(chunksize > option_copy_buffer_size) {
	    option_copy_buffer_size = chunksize;
	    do_realloc = 1;
	}
    }
#endif	/* USE_NETCDF4 */
    if(buf && do_realloc) {
	free(buf);
	buf = 0;
    }
    if(buf == 0) {		/* first time or needs to grow */
	buf = emalloc(option_copy_buffer_size);
	memset((void*)buf,0,option_copy_buffer_size);
    }

    /* initialize variable iteration */
    NC_CHECK(nc_get_iter(igrp, varid, option_copy_buffer_size, &iterp));

    start = (size_t *) emalloc((iterp->rank + 1) * sizeof(size_t));
    count = (size_t *) emalloc((iterp->rank + 1) * sizeof(size_t));
    /* nc_next_iter() initializes start and count on first call,
     * changes start and count to iterate through whole variable on
     * subsequent calls. */
    while((ntoget = nc_next_iter(iterp, start, count)) > 0) {
	NC_CHECK(nc_get_vara(igrp, varid, start, count, buf));
	NC_CHECK(nc_put_vara(ogrp, ovarid, start, count, buf));
#ifdef USE_NETCDF4
	/* we have to explicitly free values for strings and vlens */
	if(vartype == NC_STRING) {
	    NC_CHECK(nc_free_string(ntoget, (char **)buf));
	} else if(vartype > NC_STRING) { /* user-defined type */
	    nc_type vclass;
	    NC_CHECK(nc_inq_user_type(igrp, vartype, NULL, NULL, NULL, NULL, &vclass));
	    if(vclass == NC_VLEN) {
		NC_CHECK(nc_free_vlens(ntoget, (nc_vlen_t *)buf));
	    }
	}
#endif	/* USE_NETCDF4 */
    } /* end main iteration loop */
#ifdef USE_NETCDF4
    /* We're all done with this input and output variable, so if
     * either variable is chunked, free up its variable chunk cache */
    /* NC_CHECK(free_var_chunk_cache(igrp, varid)); */
    /* NC_CHECK(free_var_chunk_cache(ogrp, ovarid)); */
#endif	/* USE_NETCDF4 */
    free(start);
    free(count);
    NC_CHECK(nc_free_iter(iterp));
    return stat;
}

/* Copy data from variables in group igrp to variables in
 * corresponding group with parent ogrp, and all subgroups
 * recursively  */
static int
copy_data(int igrp, int ogrp)
{
    int stat = NC_NOERR;
    int ogid;
    int nvars;
    int varid;
#ifdef USE_NETCDF4
    int numgrps;
    int *grpids;
    int i;
#endif

    int iv;			/* variable number */
    idnode_t* vlist = NULL;	/* list for vars specified with -v option */

    /*
     * If any vars were specified with -v option, get list of
     * associated variable ids relative to this group.  Assume vars
     * specified with syntax like "grp1/grp2/varname" or
     * "/grp1/grp2/varname" if they are in groups.
     */
    vlist = newidlist();	/* list for vars specified with -v option */
    for (iv=0; iv < option_nlvars; iv++) {
        if(nc_inq_gvarid(igrp, option_lvars[iv], &varid) == NC_NOERR)
            idadd(vlist, varid);
    }
    
    /* get groupid in output corresponding to group igrp in input,
     * given parent group (or root group) ogrp in output */
    NC_CHECK(get_grpid(igrp, ogrp, &ogid));
    
    /* Copy data from this group */
    NC_CHECK(nc_inq_nvars(igrp, &nvars));

    for (varid = 0; varid < nvars; varid++) {
	if (option_nlvars > 0 && ! idmember(vlist, varid))
            continue;
        if (!group_wanted(igrp, option_nlgrps, option_grpids))
            continue;
	NC_CHECK(copy_var_data(igrp, varid, ogid));
    }
#ifdef USE_NETCDF4
    /* Copy data from subgroups */
    stat = nc_inq_grps(igrp, &numgrps, NULL);
    grpids = (int *)emalloc((numgrps + 1) * sizeof(int));
    NC_CHECK(nc_inq_grps(igrp, &numgrps, grpids));

    for(i = 0; i < numgrps; i++) {
        if (!option_grpstruct && !group_wanted(grpids[i], option_nlgrps, option_grpids))
            continue;
	NC_CHECK(copy_data(grpids[i], ogid));
    }
    free(grpids);
#endif	/* USE_NETCDF4 */
    freeidlist(vlist);
    return stat;
}

/* Count total number of dimensions in ncid and all its descendant subgroups */
int
count_dims(ncid) {
    int numgrps;
    int ndims;
    NC_CHECK(nc_inq_ndims(ncid, &ndims));
#ifdef USE_NETCDF4
    NC_CHECK(nc_inq_grps(ncid, &numgrps, NULL));
    if(numgrps > 0) {
	int igrp;
	int *grpids = emalloc(numgrps * sizeof(int));
	NC_CHECK(nc_inq_grps(ncid, &numgrps, grpids));
	for(igrp = 0; igrp < numgrps; igrp++) {
	    ndims += count_dims(grpids[igrp]);
	}
	free(grpids); 
    }
#endif	/* USE_NETCDF4 */
    return ndims;
}

/* Test if special case: netCDF-3 file with more than one record
 * variable.  Performance can be very slow for this case when the disk
 * block size is large, there are many record variables, and a
 * record's worth of data for some variables is smaller than the disk
 * block size.  In this case, copying the record variables a variable
 * at a time causes much rereading of record data, so instead we want
 * to copy data a record at a time. */
static int
nc3_special_case(int ncid, int kind) {
    if (kind == NC_FORMAT_CLASSIC ||  kind == NC_FORMAT_64BIT) {
	int recdimid = 0;
	NC_CHECK(nc_inq_unlimdim(ncid, &recdimid));
	if (recdimid != -1) {	/* we have a record dimension */
	    int nvars;
	    int varid;
	    NC_CHECK(nc_inq_nvars(ncid, &nvars));
	    for (varid = 0; varid < nvars; varid++) {
		int *dimids = 0;
		int ndims;
		NC_CHECK( nc_inq_varndims(ncid, varid, &ndims) );
		if (ndims > 0) {
		    int dimids0;
		    dimids = (int *) emalloc((ndims + 1) * sizeof(int));
		    NC_CHECK( nc_inq_vardimid(ncid, varid, dimids) );
		    dimids0 = dimids[0];
		    free(dimids);
		    if(dimids0 == recdimid) {
			return 1; /* found a record variable */
		    }
		}
	    }
	}
    }
    return 0;
}

/* Classify variables in ncid as either fixed-size variables (with no
 * unlimited dimension) or as record variables (with an unlimited
 * dimension) */
static int
classify_vars(
    int ncid,	/* netCDF ID */
    size_t *nf,	/* for returning number of fixed-size variables */
    int **fvars,	/* the array of fixed_size variable IDS, caller should free */
    size_t *nr,	/* for returning number of record variables */
    int **rvars)	/* the array of record variable IDs, caller should free */
{
    int varid;
    int nvars;
    NC_CHECK(nc_inq_nvars(ncid, &nvars));
    *nf = 0;
    *fvars = (int *) emalloc(nvars * sizeof(int));
    *nr = 0;
    *rvars = (int *) emalloc(nvars * sizeof(int));
    for (varid = 0; varid < nvars; varid++) {
	if (isrecvar(ncid, varid)) {
	    (*rvars)[*nr] = varid;
	    (*nr)++;
	} else {
	    (*fvars)[*nf] = varid;
	    (*nf)++;
	}
    }
    return NC_NOERR;
}

/* Only called for classic format or 64-bit offset format files, to speed up special case */
static int
copy_fixed_size_data(int igrp, int ogrp, size_t nfixed_vars, int *fixed_varids) {
    size_t ivar;
    /* for each fixed-size variable, copy data */
    for (ivar = 0; ivar < nfixed_vars; ivar++) {
	int varid = fixed_varids[ivar];
	NC_CHECK(copy_var_data(igrp, varid, ogrp));
    }
    if (fixed_varids)
	free(fixed_varids);
    return NC_NOERR;
}

/* copy a record's worth of data for a variable from input to output */
static int
copy_rec_var_data(int ncid, 	/* input */
		  int ogrp, 	/* output */
		  int irec, 	/* record number */
		  int varid, 	/* input variable id */
		  int ovarid, 	/* output variable id */
		  size_t *start,   /* start indices for record data */
		  size_t *count,   /* edge lengths for record data */
		  void *buf	   /* buffer large enough to hold data */
    ) 
{
    NC_CHECK(nc_get_vara(ncid, varid, start, count, buf));
    NC_CHECK(nc_put_vara(ogrp, ovarid, start, count, buf));
    return NC_NOERR;
}

/* Only called for classic format or 64-bit offset format files, to speed up special case */
static int
copy_record_data(int ncid, int ogrp, size_t nrec_vars, int *rec_varids) {
    int unlimid;
    size_t nrecs = 0;		/* how many records? */
    size_t irec;
    size_t ivar;
    void **buf;			/* space for reading in data for each variable */
    int *rec_ovarids;		/* corresponding varids in output */
    size_t **start;
    size_t **count;
    NC_CHECK(nc_inq_unlimdim(ncid, &unlimid));
    NC_CHECK(nc_inq_dimlen(ncid, unlimid, &nrecs));
    buf = (void **) emalloc(nrec_vars * sizeof(void *));
    rec_ovarids = (int *) emalloc(nrec_vars * sizeof(int));
    start = (size_t **) emalloc(nrec_vars * sizeof(size_t*));
    count = (size_t **) emalloc(nrec_vars * sizeof(size_t*));
    /* get space to hold one record's worth of data for each record variable */
    for (ivar = 0; ivar < nrec_vars; ivar++) {
	int varid;
	int ndims;
	int *dimids;
	size_t value_size;
	int dimid;
	int ii;
	size_t nvals;
	char varname[NC_MAX_NAME];
	varid = rec_varids[ivar];
	NC_CHECK(nc_inq_varndims(ncid, varid, &ndims));
	dimids = (int *) emalloc((1 + ndims) * sizeof(int));
	start[ivar] = (size_t *) emalloc(ndims * sizeof(size_t));
	count[ivar] = (size_t *) emalloc(ndims * sizeof(size_t));
	NC_CHECK(nc_inq_vardimid (ncid, varid, dimids));
	value_size = val_size(ncid, varid);
	nvals = 1;
	for(ii = 1; ii < ndims; ii++) { /* for rec size, don't include first record dimension */
	    size_t dimlen;
	    dimid = dimids[ii];
	    NC_CHECK(nc_inq_dimlen(ncid, dimid, &dimlen));
	    nvals *= dimlen;
	    start[ivar][ii] = 0;
	    count[ivar][ii] = dimlen;
	}
	start[ivar][0] = 0;	
	count[ivar][0] = 1;	/* 1 record */
	buf[ivar] = (void *) emalloc(nvals * value_size);
	NC_CHECK(nc_inq_varname(ncid, varid, varname));
	NC_CHECK(nc_inq_varid(ogrp, varname, &rec_ovarids[ivar]));
	if(dimids)
	    free(dimids);
    }

    /* for each record, copy all variable data */
    for(irec = 0; irec < nrecs; irec++) {
	for (ivar = 0; ivar < nrec_vars; ivar++) {
	    int varid, ovarid;
	    varid = rec_varids[ivar];
	    ovarid = rec_ovarids[ivar];
	    start[ivar][0] = irec;
	    NC_CHECK(copy_rec_var_data(ncid, ogrp, irec, varid, ovarid, 
				       start[ivar], count[ivar], buf[ivar]));
	}
    }
    for (ivar = 0; ivar < nrec_vars; ivar++) {
	if(start[ivar])
	    free(start[ivar]);
	if(count[ivar])
	    free(count[ivar]);
    }
    if(start)
	free(start);
    if(count)
	free(count);
    for (ivar = 0; ivar < nrec_vars; ivar++) {
	if(buf[ivar]) {
	    free(buf[ivar]);
	}
    }
    if (rec_varids)
	free(rec_varids);
    if(buf)
	free(buf);
    if(rec_ovarids)
	free(rec_ovarids);
    return NC_NOERR;
}

/* copy infile to outfile using netCDF API
 */
static int
copy(char* infile, char* outfile)
{
    int stat = NC_NOERR;
    int igrp, ogrp;
    int inkind, outkind;
    int open_mode = NC_NOWRITE;
    int create_mode = NC_CLOBBER;
    size_t ndims;

    if(option_read_diskless) {
	open_mode |= NC_DISKLESS;
    }

    NC_CHECK(nc_open(infile, open_mode, &igrp));

    NC_CHECK(nc_inq_format(igrp, &inkind));

/* option_kind specifies which netCDF format for output: 
 *   -1 -> same as input, 
 *    1 -> classic
 *    2 -> 64-bit offset
 *    3 -> netCDF-4, 
 *    4 -> netCDF-4 classic model
 *
 * However, if compression or shuffling was specified and kind was -1,
 * kind is changed to format 4 that supports compression for input of
 * type 1 or 2.  
 */
    outkind = option_kind;
    if (option_kind == SAME_AS_INPUT) {	/* default, kind not specified */
	outkind = inkind;
	/* Deduce output kind if netCDF-4 features requested */
	if (inkind == NC_FORMAT_CLASSIC || inkind == NC_FORMAT_64BIT) { 
	    if (option_deflate_level > 0 || 
		option_shuffle_vars == NC_SHUFFLE || 
		option_chunkspec) 
	    { 
		outkind = NC_FORMAT_NETCDF4_CLASSIC;
	    }
	}
    }

#ifdef USE_NETCDF4
    if(option_chunkspec) {
	/* Now that input is open, can parse option_chunkspec into binary
	 * structure. */
	NC_CHECK(chunkspec_parse(igrp, option_chunkspec));
    }
#endif	/* USE_NETCDF4 */

	/* Check if any vars in -v don't exist */
    if(missing_vars(igrp, option_nlvars, option_lvars))
	exit(EXIT_FAILURE);

    if(option_nlgrps > 0) {
	if(inkind != NC_FORMAT_NETCDF4) {
	    error("Group list (-g ...) only permitted for netCDF-4 file");
	    exit(EXIT_FAILURE);
	}
	/* Check if any grps in -g don't exist */
	if(grp_matches(igrp, option_nlgrps, option_lgrps, option_grpids) == 0)
	    exit(EXIT_FAILURE);
    }

    if(option_write_diskless)
	create_mode |= NC_WRITE | NC_DISKLESS; /* NC_WRITE persists diskless file on close */
    switch(outkind) {
    case NC_FORMAT_CLASSIC:
	/* nothing to do */
	break;
    case NC_FORMAT_64BIT:
	create_mode |= NC_64BIT_OFFSET;
	break;
#ifdef USE_NETCDF4
    case NC_FORMAT_NETCDF4:
	create_mode |= NC_NETCDF4;
	break;
    case NC_FORMAT_NETCDF4_CLASSIC:
	create_mode |= NC_NETCDF4 | NC_CLASSIC_MODEL;
	break;
#else
    case NC_FORMAT_NETCDF4:
    case NC_FORMAT_NETCDF4_CLASSIC:
	error("nccopy built with --disable-netcdf4, can't create netCDF-4 files");
	break;
#endif	/* USE_NETCDF4 */
    default:
	error("bad value (%d) for -k option\n", option_kind);
	break;
    }
    NC_CHECK(nc_create(outfile, create_mode, &ogrp));
    NC_CHECK(nc_set_fill(ogrp, NC_NOFILL, NULL));

#ifdef USE_NETCDF4
    /* Because types in one group may depend on types in a different
     * group, need to create all groups before defining types */
    if(inkind == NC_FORMAT_NETCDF4) {
	NC_CHECK(copy_groups(igrp, ogrp));
	NC_CHECK(copy_types(igrp, ogrp));
    }
#endif	/* USE_NETCDF4 */

    ndims = count_dims(igrp);
    NC_CHECK(dimmap_init(ndims));
    NC_CHECK(copy_schema(igrp, ogrp));
    NC_CHECK(nc_enddef(ogrp));

    /* For performance, special case netCDF-3 input or output file with record
     * variables, to copy a record-at-a-time instead of a
     * variable-at-a-time. */
    /* TODO: check that these special cases work with -v option */
    if(nc3_special_case(igrp, inkind)) {
	size_t nfixed_vars, nrec_vars;
	int *fixed_varids;
	int *rec_varids;
	NC_CHECK(classify_vars(igrp, &nfixed_vars, &fixed_varids, &nrec_vars, &rec_varids));
	NC_CHECK(copy_fixed_size_data(igrp, ogrp, nfixed_vars, fixed_varids));
	NC_CHECK(copy_record_data(igrp, ogrp, nrec_vars, rec_varids));
    } else if (nc3_special_case(ogrp, outkind)) {
	size_t nfixed_vars, nrec_vars;
	int *fixed_varids;
	int *rec_varids;
	/* classifies output vars, but returns input varids */
	NC_CHECK(classify_vars(ogrp, &nfixed_vars, &fixed_varids, &nrec_vars, &rec_varids));
	NC_CHECK(copy_fixed_size_data(igrp, ogrp, nfixed_vars, fixed_varids));
	NC_CHECK(copy_record_data(igrp, ogrp, nrec_vars, rec_varids));
    } else {	    
	NC_CHECK(copy_data(igrp, ogrp)); /* recursive, to handle nested groups */
    }

    NC_CHECK(nc_close(igrp));
    NC_CHECK(nc_close(ogrp));
    return stat;
}

/* 
 * For non-negative numeric string with multiplier suffix K, M, G, T,
 * or P (or lower-case equivalent), return corresponding value
 * incorporating multiplier 1000, 1000000, 1.0d9, ... 1.0d15, or -1.0
 * for error.
 */
static double
double_with_suffix(char *str) {
    double dval;
    char *suffix = 0;
    errno = 0;
    dval = strtod(str, &suffix);
    if(dval < 0 || errno != 0)
	return -1.0;
    if(*suffix) {
	switch (*suffix) {
	case 'k': case 'K':
	    dval *= 1000;
	    break;
	case 'm': case 'M':
	    dval *= 1000000;
	    break;
	case 'g': case 'G':
	    dval *= 1000000000;
	    break;
	case 't': case 'T':
	    dval *= 1.0e12;
	    break;
	case 'p': case 'P':
	    dval *= 1.0e15;
	    break;
	default:
	    dval = -1.0;	/* error, suffix multiplier must be K, M, G, or T */
	}		
    }
    return dval;
}

static void
usage(void)
{
#define USAGE   "\
  [-k n]    specify kind of netCDF format for output file, default same as input\n\
	    1 classic, 2 64-bit offset, 3 netCDF-4, 4 netCDF-4 classic model\n\
  [-d n]    set deflation compression level, default same as input (0=none 9=max)\n\
  [-s]      add shuffle option to deflation compression\n\
  [-c chunkspec] specify chunking for dimensions, e.g. \"dim1/N1,dim2/N2,...\"\n\
  [-u]      convert unlimited dimensions to fixed-size dimensions in output copy\n\
  [-w]      write whole output file from diskless netCDF on close\n\
  [-v var1,...] include data for only listed variables, but definitions for all variables\n\
  [-V var1,...] include definitions and data for only listed variables\n\
  [-g grp1,...] include data for only variables in listed groups, but all definitions\n\
  [-G grp1,...] include definitions and data only for variables in listed groups\n\
  [-m n]    set size in bytes of copy buffer, default is 5000000 bytes\n\
  [-h n]    set size in bytes of chunk_cache for chunked variables\n\
  [-e n]    set number of elements that chunk_cache can hold\n\
  [-r]      read whole input file into diskless file on open (classic or 64-bit offset format only)\n\
  infile    name of netCDF input file\n\
  outfile   name for netCDF output file\n"

    /* Don't document this flaky option until it works better */
    /* [-x]      use experimental computed estimates for variable-specific chunk caches\n\ */

    error("%s [-k n] [-d n] [-s] [-c chunkspec] [-u] [-w] [-[v|V] varlist] [-[g|G] grplist] [-m n] [-h n] [-e n] [-r] infile outfile\n%s",
	  progname, USAGE);
}

int
main(int argc, char**argv)
{
    char* inputfile = NULL;
    char* outputfile = NULL;
    int c;

/* table of formats for legal -k values */
    struct Kvalues {
	char* name;
	int kind;
    } legalkinds[] = {
	{"1", NC_FORMAT_CLASSIC},
	{"classic", NC_FORMAT_CLASSIC},
	
	/* The 64-bit offset kind (2) */
	{"2", NC_FORMAT_64BIT},
	{"64-bit-offset", NC_FORMAT_64BIT},
	{"64-bit offset", NC_FORMAT_64BIT},
	
	/* NetCDF-4 HDF5 format */
	{"3", NC_FORMAT_NETCDF4},
	{"hdf5", NC_FORMAT_NETCDF4},
	{"netCDF-4", NC_FORMAT_NETCDF4},
	{"netCDF4", NC_FORMAT_NETCDF4},
	{"enhanced", NC_FORMAT_NETCDF4},

	/* NetCDF-4 HDF5 format, but using only nc3 data model */
	{"4", NC_FORMAT_NETCDF4_CLASSIC},
	{"hdf5-nc3", NC_FORMAT_NETCDF4_CLASSIC},
	{"netCDF-4 classic model", NC_FORMAT_NETCDF4_CLASSIC},
	{"netCDF4_classic", NC_FORMAT_NETCDF4_CLASSIC},
	{"enhanced-nc3", NC_FORMAT_NETCDF4_CLASSIC},

	/* null terminate*/
	{NULL,0}
    };

    opterr = 1;
    progname = argv[0];

    if (argc <= 1)
    {
       usage();
    }

    while ((c = getopt(argc, argv, "k:d:sum:c:h:e:rwxg:G:v:V:")) != -1) {
	switch(c) {
        case 'k': /* for specifying variant of netCDF format to be generated 
                     Possible values are:
                     1 (=> classic 32 bit)
                     2 (=> classic 64 bit offsets)
                     3 (=> netCDF-4/HDF5)
                     4 (=> classic, but stored in netCDF-4/HDF5 format)
                     Also allow string versions of above
                     "classic"
                     "64-bit-offset"
                     "64-bit offset"
		     "enhanced" | "hdf5" | "netCDF-4"
                     "enhanced-nc3" | "hdf5-nc3" | "netCDF-4 classic model"
		   */
	    {
		struct Kvalues* kvalue;
		char *kind_name = (char *) emalloc(strlen(optarg)+1);
		(void)strcpy(kind_name, optarg);
	        for(kvalue=legalkinds;kvalue->name;kvalue++) {
		    if(strcmp(kind_name,kvalue->name) == 0) {
		        option_kind = kvalue->kind;
			break;
		    }
		}
		if(kvalue->name == NULL) {
		    error("invalid format: %s", kind_name);
		}
	    }
	    break;
	case 'd':		/* non-default compression level specified */
	    option_deflate_level = strtol(optarg, NULL, 10);
	    if(option_deflate_level < 0 || option_deflate_level > 9) {
		error("invalid deflation level: %d", option_deflate_level);
	    }
	    break;
	case 's':		/* shuffling, may improve compression */
	    option_shuffle_vars = NC_SHUFFLE;
	    break;
	case 'u':		/* convert unlimited dimensions to fixed size */
	    option_fix_unlimdims = 1;
	    break;
	case 'm':		/* non-default size of data copy buffer */
	{
	    double dval = double_with_suffix(optarg);	/* "K" for kilobytes. "M" for megabytes, ... */
	    if(dval < 0)
		error("Suffix used for '-m' option value must be K, M, G, T, or P");
	    option_copy_buffer_size = dval;
	    break;
	}
	case 'h':		/* non-default size of chunk cache */
	{
	    double dval = double_with_suffix(optarg);	/* "K" for kilobytes. "M" for megabytes, ... */
	    if(dval < 0)
		error("Suffix used for '-h' option value must be K, M, G, T, or P");
	    option_chunk_cache_size = dval;
	    break;
	}
	case 'e':		/* number of elements chunk cache can hold */
	{
	    double dval = double_with_suffix(optarg);	/* "K" for kilobytes. "M" for megabytes, ... */
	    if(dval < 0 )
		error("Suffix used for '-e' option value must be K, M, G, T, or P");
	    option_chunk_cache_nelems = (long)dval;
	    break;
	}
	case 'r':
	    option_read_diskless = 1; /* read into memory on open */
	    break;
	case 'w':
	    option_write_diskless = 1; /* write to memory, persist on close */
	    break;
	case 'x':		/* use experimental variable-specific chunk caches */
	    option_compute_chunkcaches = 1;
	    break;
	case 'c':               /* optional chunking spec for each dimension in list */
	    /* save chunkspec string for parsing later, once we know input ncid */
	    option_chunkspec = strdup(optarg);
	    break;
	case 'g':		/* group names */
	    /* make list of names of groups specified */
	    make_lgrps (optarg, &option_nlgrps, &option_lgrps, &option_grpids);
	    option_grpstruct = true;
	    break;
	case 'G':		/* group names */
	    /* make list of names of groups specified */
	    make_lgrps (optarg, &option_nlgrps, &option_lgrps, &option_grpids);
	    option_grpstruct = false;
	    break;
	case 'v':		/* variable names */
	    /* make list of names of variables specified */
	    make_lvars (optarg, &option_nlvars, &option_lvars);
	    option_varstruct = true;
	    break;
	case 'V':		/* variable names */
	    /* make list of names of variables specified */
	    make_lvars (optarg, &option_nlvars, &option_lvars);
	    option_varstruct = false;
	    break;
	default: 
	    usage();
        }
    }
    argc -= optind;
    argv += optind;

    if (argc != 2) {
	error("one input file and one output file required");
    }
    inputfile = argv[0];
    outputfile = argv[1];

    if(strcmp(inputfile, outputfile) == 0) {
	error("output would overwrite input");
    }

    if(copy(inputfile, outputfile) != NC_NOERR)
        exit(EXIT_FAILURE);
    exit(EXIT_SUCCESS);
}
END_OF_MAIN();
