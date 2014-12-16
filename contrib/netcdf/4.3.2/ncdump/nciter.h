/*********************************************************************
 *   Copyright 2009, University Corporation for Atmospheric Research
 *   See netcdf/README file for copying and redistribution conditions.
 *   "$Id: nciter.h 400 2010-08-27 21:02:52Z russ $"
 *********************************************************************/

#ifndef _NCITER_
#define _NCITER_

#include <netcdf.h>

#if defined(__cplusplus)
extern "C" {
#endif

/* 
 * The opaque structure to hold per-variable state of data iteration
 */
typedef struct {
    int first;	     /* false after first invocation of nc_next_iter() */
    int right_dim;   /* rightmost dimension for start of variable pieces */
    size_t rows;     /* how many subpiece rows in bufsiz */
    size_t numrows;  /* how many row pieces in right_dim dimension */
    size_t cur;	     /* current "row" in loop over row pieces */
    size_t leftover; /* how many rows left over after partitioning
		      * bufsiz into subpiece blocks */
    int more;	     /* whether there is more data still to get */
    size_t to_get;   /* number of values to get on this access */
    int rank;	     /* number of dimensions */
    size_t inc;	     /* increment for right_dim element of start vector */
    int chunked;     /* 1 if chunked, 0 if contiguous */
    size_t *dimsizes;
    size_t *chunksizes; /* ignored if not chunked */
} nciter_t;

/*
 * The Interfaces
 */

/* Get iterator for variable data.  Returns pointer to malloc'd
 * nciter_t, which caller must later release using nc_free_iter(), not
 * free(). */
extern int
nc_get_iter(int ncid, int varid, size_t bufsize, nciter_t **iterpp);

/* Iterate over blocks of variable values, using start and count
 * vectors.  Returns number of values to access (product of counts),
 * or 0 if done. */
extern size_t
nc_next_iter(nciter_t *iterp, size_t *start, size_t *count);

/* Release memory allocated for iterator */
extern int
nc_free_iter(nciter_t *iterp);

#if defined(__cplusplus)
}
#endif

#endif /* _NCITER_ */
