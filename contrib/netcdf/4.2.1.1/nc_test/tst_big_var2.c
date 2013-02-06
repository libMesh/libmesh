/*
  Copyright 2008, UCAR/Unidata
  See COPYRIGHT file for copying and redistribution conditions.

  This program tests the fix for a large file bug in versions previous
  to netCDF-4.1.2 for 32-bit platforms, writing to a variable with
  more than 1 dimension and more than 2**32 values, where the write
  starts after the first 2**32 elements.

  $Id: tst_big_var2.c,v 1.3 2010/05/19 16:38:44 russ Exp $
*/

#include <nc_tests.h>
#include <netcdf.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define FILE_NAME "tst_big_var2.nc"

/* Test with both classic and 64-bit offset files. If netcdf-4 is
 * included, test with both netCDF-4 format variants also. */
#ifdef USE_NETCDF4
#define NUM_FORMATS (4)
#else
#define NUM_FORMATS (2)
#endif

#define NUMDIMS 3		/* rank of variable in tests */
#define DIM0 2149		/* just big enough to demonstrate bug */
#define DIM1 1000
#define DIM2 2000		/* DIM0*DIM1*DIM2 > 2**32 */

/* 
 * This program tests the fix for a large file bug in versions
 * previous to netCDF-4.1.2 for 32-bit platforms, writing to a
 * variable with more than 1 dimension and more than 2**32 values,
 * where the write starts after the first 2**32 elements.  The bug
 * applies to record variables with more than 2**32 values per record
 * as well, but that's not tested here.
 */
static int
test_big_var(const char *testfile) 
{
    int ncid, varid, dimids[NUMDIMS];
    size_t index[NUMDIMS];
    int nval = 99;
    int nval_in;
    size_t start[NUMDIMS] = {0, 0, 0};
    size_t count[NUMDIMS] = {1, DIM1, DIM2};
    signed char data[DIM1][DIM2];
    int i, j;
    int nerrs = 0;
    
    /* Create a file with one big variable. */
    if (nc_create(testfile, NC_CLOBBER, &ncid)) ERR;
    if (nc_set_fill(ncid, NC_NOFILL, NULL)) ERR;
    if (nc_def_dim(ncid, "dim0", DIM0, &dimids[0])) ERR;
    if (nc_def_dim(ncid, "dim1", DIM1, &dimids[1])) ERR;
    /* if (nc_def_dim(ncid, "dim2", DIM2 - 1, &dimids[1])) ERR; */
    if (nc_def_dim(ncid, "dim2", DIM2, &dimids[2])) ERR;
    if (nc_def_var(ncid, "var", NC_BYTE, NUMDIMS, dimids, &varid)) ERR;
    if (nc_enddef(ncid)) ERR;

    /* Initialize slab of data. */
    for (i = 0; i < DIM1; i++)
	for (j = 0; j < DIM2; j++) {
	    data[i][j] = 42;
	}
    /* Just write the first and last slabs */
    start[0] = 0;
    if (nc_put_vara_schar(ncid, varid, start, count, &data[0][0])) ERR;
    for (i = 0; i < DIM1; i++)
	for (j = 0; j < DIM2; j++) {
	    data[i][j] = 19;
	}
    start[0] = DIM0 - 1;
    if (nc_put_vara_schar(ncid, varid, start, count, &data[0][0])) ERR;
    if (nc_close(ncid)) ERR;

    /* Open the file and check it. */
    if (nc_open(testfile, NC_NOWRITE, &ncid)) ERR;
    if (nc_inq_varid(ncid, "var", &varid)) ERR;
    /* Read and check data in the first and last slabs */
    start[0] = 0;
    if (nc_get_vara_schar(ncid, varid, start, count, &data[0][0])) ERR;
    for (i = 0; i < DIM1; i++)
	for (j = 0; j < DIM2; j++)
	{
	    if (data[i][j] != 42 )
	    {
		printf("error on start[0]: %d i: %d j: %d expected %d got %d\n", 
		       start[0], i, j, 42, data[i][j]);
		ERR;
		if(nerrs++ > 1)
		    return;
	    }
	}
    start[0] = DIM0 - 1;
    if (nc_get_vara_schar(ncid, varid, start, count, &data[0][0])) ERR;
    for (i = 0; i < DIM1; i++)
	for (j = 0; j < DIM2; j++)
	{
	    if (data[i][j] != 19 )
	    {
		printf("error on start[0]: %d i: %d j: %d expected %d got %d\n", 
		       start[0], i, j, 19, data[i][j]);
		ERR;
		if(nerrs++ > 1)
		    return;
	    }
	}
    if (nc_close(ncid)) ERR;
    return 0;
}

int
main(int argc, char **argv) {
    int i;
    char testfile[NC_MAX_NAME + 1];

    printf("\n*** Testing multidimensional variable with more than 2**32 values\n");
    sprintf(testfile, "%s/%s", TEMP_LARGE, FILE_NAME);
    for (i = NC_FORMAT_CLASSIC; i <= NUM_FORMATS; i++)
    {
       printf("*** testing format %d file with byte variable with > 2**32 values...", i);
       nc_set_default_format(i, NULL);
       test_big_var(testfile);
       (void) remove(testfile);
       SUMMARIZE_ERR;
   }
    
    FINAL_RESULTS;
}
