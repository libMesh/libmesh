/* This is part of the netCDF package.
   Copyright 2005 University Corporation for Atmospheric Research/Unidata
   See COPYRIGHT file for conditions of use.
*/

/*
 * Example illustrates the use of SZIP compression in netCDF5
 * Taken from HDF5 example.
 */

#include <config.h>
#include <stdlib.h>
#include <stdio.h>
#include "nc_tests.h"
#include "err_macros.h"
#include "netcdf.h"

#undef PLAIN
#define USECLOSE

/* Szip Constants. */
#define HDF5_FILTER_SZIP 4
#define H5_SZIP_MAX_PIXELS_PER_BLOCK_IN 32
#define H5_SZIP_MAX_PIXELS_PER_BLOCK_OUT 64

/* Option Mask Flags (Refere to HDF5 szip documentation)  */
#define H5_SZIP_ALLOW_K13_OPTION_MASK 1 /*Allows k split = 13 compression mode. (Default)*/
#define H5_SZIP_CHIP_OPTION_MASK      2 /*Compresses exactly as in hardware*/
#define H5_SZIP_EC_OPTION_MASK        4 /*Selects entropy coding method. (Default)*/
#define H5_SZIP_NN_OPTION_MASK       32 /*Selects nearest neighbor coding method*/

#define NX 500
#define NY 600
#define CH_NX 100
#define CH_NY 25


static void initialize(void);
static int compare(void);

static float buf[NX][NY];
static float buf_r[NX][NY];

int
main(void)
{
    int ncid, varid, dimids[2];
    size_t dims[2], chunk_size[2];
    unsigned int szip_params[2]; /* [0]=options_mask [1]=pixels_per_block */
    int errcnt = 0;

    /* Create a new file using read/write access. */
    if(nc_create("testszip.nc", NC_CLOBBER|NC_NETCDF4, &ncid)) ERR;

    /* Create dims */
    dims[0] = NX;
    dims[1] = NY;
    if(nc_def_dim(ncid, "x", dims[0], &dimids[0])) ERR;
    if(nc_def_dim(ncid, "y", dims[1], &dimids[1])) ERR;

    /* Create a dimensioned variable */
    if(nc_def_var(ncid, "datasetF32", NC_FLOAT, 2, dimids, &varid)) ERR;

    /* no fill */
    if(nc_def_var_fill(ncid, varid, 1, NULL)) ERR;

    /* Define chunking for the variable:
    * the raw data is to be partitioned into 100x100 element chunks.
    */
    chunk_size[0] = CH_NX;
    chunk_size[1] = CH_NY;
    if(nc_def_var_chunking(ncid, varid, NC_CHUNKED, chunk_size)) ERR;

#ifndef PLAIN
    /*
     * Set parameters for SZIP compression; check the description of
     * the H5Pset_szip function in the HDF5 Reference Manual for more
     * information.
     */
    szip_params[0] = H5_SZIP_NN_OPTION_MASK;
    szip_params[1] = H5_SZIP_MAX_PIXELS_PER_BLOCK_IN;
    { int stat = nc_def_var_filter(ncid, varid, HDF5_FILTER_SZIP, 2, szip_params);
       if(stat) {
         fprintf(stderr,"XXX: %d %s\b",stat,nc_strerror(stat));
         ERR;
       }
    }
#endif

    if(nc_enddef(ncid)) ERR;

    initialize();

    /* Write the array to the file */
    if(nc_put_var_float(ncid, varid, &buf[0][0])) ERR;

#ifdef USECLOSE
    /* Close and re-open the file */
    if(nc_close(ncid)) ERR;
    if(nc_open("testszip.nc", NC_NETCDF4, &ncid)) ERR;
    if(nc_inq_varid(ncid, "datasetF32", &varid)) ERR;
#endif

    /*
    * Read the array.  This is similar to writing data,
    * except the data flows in the opposite direction.
    * Note: Decompression should be automatic.
    */

    memset(buf_r,0,sizeof(buf_r));
    if(nc_get_var_float(ncid, varid, &buf_r[0][0])) ERR;

    /* Do comparison */
    errcnt = compare();

    if(nc_close(ncid)) ERR;

    if(errcnt) ERR;

    SUMMARIZE_ERR;
    FINAL_RESULTS;

    return (errcnt==0?0:1);
}

static int
compare()
{
    int i,j;
    int errs = 0;
    
    /* Do comparison */
    for (i=0; i < NX; i++) {
	for (j=0; j < NY; j++) {
            if(buf[i][j] != buf_r[i][j]) {
		errs++;
	        printf("mismatch: [%d][%d]: write = %f read=%f\n",
		        i,j,buf[i][j],buf_r[i][j]);
	}
     }
   }
   return errs;
}

static void
initialize(void)
{
   int i, j;
   /* Initialize data buffer with some bogus data. */
   for(i=0; i < NX; i++) {
     for(j=0; j < NY; j++) {
       buf[i][j] = (float)(i + j);
     }
   }
}

