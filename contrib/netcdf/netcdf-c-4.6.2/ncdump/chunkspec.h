/*********************************************************************
 *   Copyright 2010, UCAR/Unidata
 *   See netcdf/COPYRIGHT file for copying and redistribution conditions.
 *   $Id $
 *********************************************************************/
#ifndef _CHUNKSPEC_H_
#define _CHUNKSPEC_H_

/* Parse chunkspec string and convert into internal data structure,
 * associating dimids from open file or group specified by ncid with
 * corresponding chunk sizes */
extern int
chunkspec_parse(int ncid, const char *spec);

/* Return chunk size in chunkspec string specified for dimension
 * corresponding to dimid, 0 if not found */
extern size_t
dimchunkspec_size(int dimid);

/* Return number of dimensions for which chunking was specified in
 * chunkspec string on command line, 0 if no chunkspec string was
 * specified. */
extern int
dimchunkspec_ndims(void);

/* Return whether chunking should be omitted, due to explicit
 * command-line specification. */
extern bool_t
dimchunkspec_omit(void);

extern bool_t varchunkspec_omit(int grpid, int varid);

extern size_t* varchunkspec_chunksizes(int grpid, int varid);

extern size_t varchunkspec_ndims(int grpid, int varid);

extern bool_t varchunkspec_exists(int grpid, int varid);

extern void chunkspecinit(void);


#endif	/* _CHUNKSPEC_H_  */
