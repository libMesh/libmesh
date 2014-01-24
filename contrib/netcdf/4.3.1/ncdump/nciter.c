/*********************************************************************
 *   Copyright 2009, University Corporation for Atmospheric Research
 *   See netcdf/README file for copying and redistribution conditions.
 *   "$Id: nciter.c 400 2010-08-27 21:02:52Z russ $"
 *********************************************************************/

#include "config.h"		/* for USE_NETCDF4 macro */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <netcdf.h>
#include "utils.h"
#include "nciter.h"

/* Initialize block iteration for variables, including those that
 * won't fit in the copy buffer all at once.  Returns error if
 * variable is chunked but size of chunks is too big to fit in bufsize
 * bytes. */
static int
nc_blkio_init(size_t bufsize, 	/* size in bytes of in-memory copy buffer */
	      size_t value_size, /* size in bytes of each variable element */
	      int rank,		 /* number of dimensions for variable */
	      int chunked,	  /* 1 if variable is chunked, 0 otherwise */
	      nciter_t *iter /* returned iteration state, don't mess with it */
    ) {
    int stat = NC_NOERR;
    int i;
    long long prod;
    size_t *dims = iter->dimsizes;

    iter->rank = rank;
    iter->first = 1;
    iter->more = 1;
    iter->chunked = chunked;
    prod = value_size;
    if(iter->chunked == 0) {	/* contiguous */
	iter->right_dim = rank - 1;
	for(i = rank; i > 0; i--) {
	    if(prod*dims[i-1] <= bufsize) {
		prod *= dims[i-1];
		iter->right_dim--;
	    } else {
		break;
	    }
	}
	if (i > 0) {	     /* variable won't fit in bufsize bytes */
	    iter->rows = bufsize/prod;
	    iter->numrows = dims[iter->right_dim] / iter->rows;
	    iter->leftover = dims[iter->right_dim] - iter->numrows * iter->rows;
	    iter->cur = 1;
	    iter->inc = iter->rows;
	    return stat;
	}
	/* else, variable will fit in bufsize bytes of memory. */
	iter->right_dim = 0;
	iter->rows = dims[0];
	iter->inc = 0;
	return stat;
    } 
    /* else, handle chunked case */
    for(i = 0; i < rank; i++) {
	prod *= iter->chunksizes[i];
    }
    if(prod > bufsize) {
	stat = NC_ENOMEM;
	fprintf(stderr, "chunksize (= %ld) > copy_buffer size (= %ld)\n",
			(long)prod, (long)bufsize);
    }
    return stat;
}

/* From netCDF type in group igrp, get size in memory needed for each
 * value.  Wouldn't be needed if nc_inq_type() was a netCDF-3 function
 * too. */
static int
inq_value_size(int igrp, nc_type vartype, size_t *value_sizep) {
    int stat = NC_NOERR;
    
#ifdef USE_NETCDF4
    NC_CHECK(nc_inq_type(igrp, vartype, NULL, value_sizep));
#else
    switch(vartype) {
    case NC_BYTE:
	*value_sizep = sizeof(signed char);
	break;
    case NC_CHAR:
	*value_sizep = sizeof(char);
	break;
    case NC_SHORT:
	*value_sizep = sizeof(short);
	break;
    case NC_INT:
	*value_sizep = sizeof(int);
	break;
    case NC_FLOAT:
	*value_sizep = sizeof(float);
	break;
    case NC_DOUBLE:
	*value_sizep = sizeof(double);
	break;
    default:
	NC_CHECK(NC_EBADTYPE);
	break;
    }
#endif	/* USE_NETCDF4 */
    return stat;
}

/*
 * Updates a vector of size_t, odometer style.  Returns 0 if odometer
 * overflowed, else 1.
 */
static int
up_start(
     int ndims,		 /* Number of dimensions */
     const size_t *dims, /* The "odometer" limits for each dimension */
     int incdim,	 /* the odmometer increment dimension */
     size_t inc,	 /* the odometer increment for that dimension */
     size_t* odom	 /* The "odometer" vector to be updated */
     )
{
    int id;
    int ret = 1;

    if(inc == 0) {
	return 0;
    }
    odom[incdim] += inc;
    for (id = incdim; id > 0; id--) {
	if(odom[id] >= dims[id]) {
	    odom[id-1]++;
	    odom[id] -= dims[id];
	}
    }
    if (odom[0] >= dims[0])
      ret = 0;
    return ret;
}

/*
 * Updates a vector of size_t, odometer style, for chunk access.
 * Returns 0 if odometer overflowed, else 1.
 */
static int
up_start_by_chunks(
     int ndims,		 /* Number of dimensions */
     const size_t *dims, /* The "odometer" limits for each dimension */
     const size_t *chunks, /* the odometer increments for each dimension */
     size_t* odom	 /* The "odometer" vector to be updated */
     )
{
    int incdim = ndims - 1;
    int id;
    int ret = 1;

    odom[incdim] += chunks[incdim];
    for (id = incdim; id > 0; id--) {
	if(odom[id] >= dims[id]) {
	    odom[id-1] += chunks[id-1];
	    /* odom[id] -= dims[id]; */
	    odom[id] = 0;
	}
    }
    if (odom[0] >= dims[0])
      ret = 0;
    return ret;
}

/* Begin public interfaces */

/* Initialize iteration for a variable.  Just a wrapper for
 * nc_blkio_init() that makes the netCDF calls needed to initialize
 * lower-level iterator. */
int
nc_get_iter(int ncid,
	     int varid,
	     size_t bufsize,   /* size in bytes of memory buffer */
	     nciter_t **iterpp /* returned opaque iteration state */) 
{
    int stat = NC_NOERR;
    nciter_t *iterp;
    nc_type vartype;
    size_t value_size;      /* size in bytes of each variable element */
    int ndims;		    /* number of dimensions for variable */
    int *dimids;
    long long nvalues = 1;
    int dim;
    int chunked = 0;

    /* Caller should free this by calling nc_free_iter(iterp) */
    iterp = (nciter_t *) emalloc(sizeof(nciter_t));
    memset((void*)iterp,0,sizeof(nciter_t)); /* make sure it is initialized */

    NC_CHECK(nc_inq_varndims(ncid, varid, &ndims));
    
    dimids = (int *) emalloc((ndims + 1) * sizeof(int));
    
    iterp->dimsizes = (size_t *) emalloc((ndims + 1) * sizeof(size_t));
    iterp->chunksizes = (size_t *) emalloc((ndims + 1) * sizeof(size_t));

    NC_CHECK(nc_inq_vardimid (ncid, varid, dimids));
    for(dim = 0; dim < ndims; dim++) {
	size_t len;
	NC_CHECK(nc_inq_dimlen(ncid, dimids[dim], &len));
	nvalues *= len;
	iterp->dimsizes[dim] = len;
    }
    NC_CHECK(nc_inq_vartype(ncid, varid, &vartype));
    NC_CHECK(inq_value_size(ncid, vartype, &value_size));
#ifdef USE_NETCDF4    
    {
	int contig = 1;
	if(ndims > 0) {
	    NC_CHECK(nc_inq_var_chunking(ncid, varid, &contig, NULL));
	}
	if(contig == 0) {	/* chunked */
	    NC_CHECK(nc_inq_var_chunking(ncid, varid, &contig, iterp->chunksizes));
	    chunked = 1;
	}
    }
#endif	/* USE_NETCDF4 */
    NC_CHECK(nc_blkio_init(bufsize, value_size, ndims, chunked, iterp));
    iterp->to_get = 0;
    free(dimids);
    *iterpp = iterp;
    return stat;
}

/* Iterate on blocks for variables, by updating start and count vector
 * for next vara call.  Assumes nc_get_iter called first.  Returns
 * number of variable values to get, 0 if done, negative number if
 * error, so use like this: 
   size_t to_get;
   while((to_get = nc_next_iter(&iter, start, count)) > 0) { 
      ... iteration ... 
   } 
   if(to_get < 0) { ... handle error ... }
 */
size_t
nc_next_iter(nciter_t *iter,	/* returned opaque iteration state */
	     size_t *start, 	/* returned start vector for next vara call */
	     size_t *count	/* returned count vector for next vara call */
    ) {
    int i;
    /* Note: special case for chunked variables is just an
     * optimization, the contiguous code below is OK even
     * for chunked variables, but in general will do more I/O ... */
    if(iter->first) {
	if(!iter->chunked) { 	/* contiguous storage */
	    for(i = 0; i < iter->right_dim; i++) {
		start[i] = 0;
		count[i] = 1;
	    }
	    start[iter->right_dim] = 0;
	    count[iter->right_dim] = iter->rows;
	    for(i = iter->right_dim + 1; i < iter->rank; i++) {
		start[i] = 0;
		count[i] = iter->dimsizes[i];
	    }
	} else {		/* chunked storage */
	    for(i = 0; i < iter->rank; i++) {
		start[i] = 0;
		if(iter->chunksizes[i] <= iter->dimsizes[i])
		    count[i] = iter->chunksizes[i];
		else /* can happen for variables with only unlimited dimensions */
		    count[i] = iter->dimsizes[i];
	    }
	}
	iter->first = 0;
    } else {
	if(!iter->chunked) { 	/* contiguous storage */
	    iter->more = up_start(iter->rank, iter->dimsizes, iter->right_dim, 
				  iter->inc, start);
	    /* iterate on pieces of variable */
	    if(iter->cur < iter->numrows) {
		iter->inc = iter->rows;
		count[iter->right_dim] = iter->rows;
		iter->cur++;
	    } else {
		if(iter->leftover > 0) {
		    count[iter->right_dim] = iter->leftover;
		    iter->inc = iter->leftover;
		    iter->cur = 0;
		}
	    }
	} else {		/* chunked storage */
	    iter->more = up_start_by_chunks(iter->rank, iter->dimsizes, 
					    iter->chunksizes, start);
	    /* adjust count to stay in range of dimsizes */
	    for(i = 0; i < iter->rank; i++) {
		int leftover = iter->dimsizes[i] - start[i];
		if(iter->chunksizes[i] <= iter->dimsizes[i])
		    count[i] = iter->chunksizes[i];
		else /* can happen for variables with only unlimited dimensions */
		    count[i] = iter->dimsizes[i];
		if(leftover < count[i]) 
		    count[i] = leftover;
	    }
	}
    }
    iter->to_get = 1;
    for(i = 0; i < iter->rank; i++) {
	iter->to_get *= count[i];
    }
    return iter->more == 0 ? 0 : iter->to_get ;
}

/* Free iterator and its internally allocated memory */
int
nc_free_iter(nciter_t *iterp) {
    if(iterp->dimsizes)
	free(iterp->dimsizes);
    if(iterp->chunksizes)
	free(iterp->chunksizes);
    free(iterp);
    return NC_NOERR;
}
