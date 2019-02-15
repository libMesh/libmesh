/*
  Copyright 2018, UCAR/Unidata
  See COPYRIGHT file for copying and redistribution conditions.
*/
/*
This file is the same as nc_test4/test_filter.c 
*/

/*! \file
Example program for write then read of a variable using bzip2 compression.

\ingroup tutorial

This is an example which 
creates a file with a variable that is compressed using bzip2.
Then it reads that file and verifies that it returned the correct
uncompressed data.

The meta-data (.cdl) for the created file is as follows:
\code
netcdf bzip2 {
dimensions:
	dim0 = 4 ;
	dim1 = 4 ;
	dim2 = 4 ;
	dim3 = 4 ;
variables:
	float var(dim0, dim1, dim2, dim3) ;
		var:_Storage = "chunked" ;
		var:_ChunkSizes = 4, 4, 4, 4 ;
		var:_Filter = "307,9" ;
		var:_NoFill = "true" ;
data:

 var =
  0, 1, 2, 3,
  4, 5, 6, 7,
  ...
  252, 253, 254, 255 ;
}
\endcode
*/

#include "config.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <hdf5.h>
#include "netcdf.h"

/* The HDF assigned id for bzip compression */
#define BZIP2_ID 307
/* The compression level used in this example */
#define BZIP2_LEVEL 9

#define TESTFILE "bzip2.nc"

/* Point at which we give up */
#define MAXERRS 8

#define NDIMS 4
#define DIMSIZE 4
#define CHUNKSIZE 4 /* Note: not the total size of the chunk, but size wrt a dim*/

static size_t dimsize = DIMSIZE;
static size_t chunksize = CHUNKSIZE;
static size_t actualdims = NDIMS;

static size_t actualproduct = 1; /* x-product over dim sizes */
static size_t chunkproduct = 1; /* x-product over chunksizes */

static size_t dims[NDIMS];
static size_t chunks[NDIMS];

static int nerrs = 0;

static int ncid, varid;
static int dimids[NDIMS];
static float* array = NULL;
static float* expected = NULL;
static unsigned int filterid = 0;
static unsigned int* params = NULL;

/* Forward */
static void init(int argc, char** argv);
static int test_bzip2(void);
static int verifychunks(void);

#define ERRR do { \
fflush(stdout); /* Make sure our stdout is synced with stderr. */ \
fprintf(stderr, "Sorry! Unexpected result, %s, line: %d\n", \
	__FILE__, __LINE__);				    \
nerrs++;\
} while (0)

static int
check(int err,int line)
{
    if(err != NC_NOERR) {
	fprintf(stderr,"fail (%d): %s\n",line,nc_strerror(err));
	fflush(stderr);
	exit(1);
    }
    return NC_NOERR;
}

#define CHECK(x) check(x,__LINE__)

/*
Read the chunking information about the variable
and verify that it is as expected.
*/

static int
verifychunks(void)
{
    int i;
    int store = -1;
    size_t chunksizes[NDIMS];
    memset(chunksizes,0,sizeof(chunksizes));
    CHECK(nc_inq_var_chunking(ncid, varid, &store, chunksizes));
    /* Storate must be chunked, not contiguous */
    if(store != NC_CHUNKED) {
	fprintf(stderr,"bad chunk store\n");
	return NC_ESTORAGE;
    }
    /* Chunk sizes must match our predefined set */
    for(i=0;i<actualdims;i++) {
        if(chunksizes[i] != chunks[i]) {
	    fprintf(stderr,"bad chunk size: %d\n",i);
	    return NC_EBADCHUNK;
	}
    }
    return 1;
}

/*
Compare the data we wrote against the data we read.
*/

static int
compare(void)
{
    int errs = 0;
    int i;
    printf("data comparison: |array|=%ld\n",(unsigned long)actualproduct);
    for(i=0;i<actualproduct;i++) {
	if(expected[i] != array[i]) {
	    printf("mismatch: array[%d]=%f expected[%d]=%f\n",
                            i,array[i],i,expected[i]);
            errs++;
            if(errs >= MAXERRS)
                break;
	}
   }
   if(errs == 0)
        printf("no data errors\n");
   if(actualproduct <= 1)
	return NC_EBADDIM;
   return (errs == 0 ? NC_NOERR: NC_EINVAL);
}

/*
Create the file, write it, then re-read for comparison.
*/
static int
test_bzip2(void)
{
    int i;
    unsigned int level = BZIP2_LEVEL;
    unsigned int id=0;
    size_t nparams = 0;

    printf("\n*** Testing API: bzip2 compression.\n");

    /* Clear the data array */
    memset(array,0,sizeof(float)*actualproduct);

    /* Create a file */
    CHECK(nc_create(TESTFILE, NC_NETCDF4|NC_CLOBBER, &ncid));

    /* Do not use fill for this file */
    CHECK(nc_set_fill(ncid, NC_NOFILL, NULL));

    /* Define the dimensions */
    for(i=0;i<actualdims;i++) {
	char dimname[1024];
	snprintf(dimname,sizeof(dimname),"dim%d",i);
        CHECK(nc_def_dim(ncid, dimname, dims[i], &dimids[i]));
    }

    /* Define the variable */
    CHECK(nc_def_var(ncid, "var", NC_FLOAT, actualdims, dimids, &varid));

    /* Set chunking on the variable */
    CHECK(nc_def_var_chunking(ncid,varid,NC_CHUNKED,chunks));

    /* Verify that chunking succeeded */
    if(!verifychunks())
	return NC_EINVAL;
    /* Set bzip2 compression for the variable: takes one parameter == level */
    CHECK(nc_def_var_filter(ncid,varid,BZIP2_ID,1,&level));

    /* Read back the compression info and verify it */
    level = 0;
    CHECK(nc_inq_var_filter(ncid,varid,&id,&nparams,&level));
    if(id != BZIP2_ID || nparams != 1 || level != BZIP2_LEVEL) {
        printf("test_filter: filter def/inq mismatch\n");
	return NC_EFILTER;
    }
    /* Show the level */
    printf("show parameters for bzip2: level=%u\n",level);
    /* Show chunking */ 
    printf("show chunks:");
    for(i=0;i<actualdims;i++)
	printf("%s%ld",(i==0?" chunks=":","),(unsigned long)chunks[i]);
    printf("\n");

    /* prepare to write */
    CHECK(nc_enddef(ncid));

    /* Fill in the array */
    for(i=0;i<actualproduct;i++)
	expected[i] = (float)i;

    /* write array */
    CHECK(nc_put_var(ncid,varid,expected));

    /* Close file */
    CHECK(nc_close(ncid));

    /* Now re-open and verify */
    printf("\n*** Testing API: bzip2 decompression.\n");

    /* Clear the data array */
    memset(array,0,sizeof(float)*actualproduct);

    /* Open the file */
    CHECK(nc_open(TESTFILE, NC_NOWRITE, &ncid));

    /* Get the variable id */
    CHECK(nc_inq_varid(ncid, "var", &varid));

    /* Check the compression algorithm */
    filterid = 0;
    nparams = 0;
    params = NULL;
    CHECK(nc_inq_var_filter(ncid,varid,&filterid,&nparams,NULL));
    if(nparams > 0) {
        params = (unsigned int*)malloc(sizeof(unsigned int)*nparams);
	if(params == NULL)
	    return NC_ENOMEM;
        CHECK(nc_inq_var_filter(ncid,varid,&filterid,&nparams,params));
    }
    if(filterid != BZIP2_ID) {
         printf("Bzip2 id mismatch: %d\n",filterid);
	return NC_EFILTER;
    }
    if(nparams != 1 && params != NULL && params[0] != BZIP2_LEVEL) {
	printf("Compression parameter mismatch\n");
	return NC_EFILTER; 
    }

    /* Verify chunking */
    if(!verifychunks())
	return 0;

    /* Read the data */
    CHECK(nc_get_var_float(ncid, varid, array));

    /* Close the file */
    CHECK(nc_close(ncid));
    return (compare() == NC_NOERR ? 0 : 1);
}

/**************************************************/
/* Utilities */

static void
init(int argc, char** argv)
{
    int i;
    /* Setup various variables */
    actualproduct = 1;
    chunkproduct = 1;
    for(i=0;i<NDIMS;i++) {
	dims[i] = dimsize;
	chunks[i] = chunksize;
	if(i < actualdims) {
	    actualproduct *= dims[i];
	    chunkproduct *= chunks[i];
	}
    }
    /* Allocate max size */
    array = (float*)calloc(1,sizeof(float)*actualproduct);
    expected = (float*)calloc(1,sizeof(float)*actualproduct);
}

/**************************************************/
int
main(int argc, char **argv)
{
    H5Eprint(stderr);
    init(argc,argv);
    if(test_bzip2() != NC_NOERR) ERRR;
    exit(nerrs > 0?1:0);
}

