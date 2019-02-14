/* This is part of the netCDF package. Copyright 2010 University
   Corporation for Atmospheric Research/Unidata.  See COPYRIGHT file for
   conditions of use. See www.unidata.ucar.edu for more info.

   Create a chunkable test file for nccopy to test chunking.
*/

#include <netcdf.h>
#include <nc_tests.h>
#include "err_macros.h"

#define DEBUG

static int ret = NC_NOERR;

/* Make trackable ERR macro replacement */
static int lerr(int stat, const char* file, int lineno) {
    fflush(stdout); /* Make sure our stdout is synced with stderr. */
    err++;
    fprintf(stderr, "Sorry! Unexpected result(%d), %s, line: %d\n",ret,file,lineno);
    fflush(stderr);                                             \
    return 2;                                                   \
}
#define LERR lerr(ret,__FILE__,__LINE__)

#define FILE_NAME "tst_chunking.nc"
#define VAR_RANK 7
#define IVAR_NAME "ivar"
#define FVAR_NAME "fvar"
#define GRP_NAME "g"
#define UNLIM_NAME "unlimited"
#define UNLIM_SIZE 10
#define DEFLATE_LEVEL 1
#define NVALS 45360		/* 7 * 4 * 2 * 3 * 5 * 6 * 9 */

static const char *dim_names[VAR_RANK] = {"dim0", "dim1", "dim2", "dim3", "dim4", "dim5", "dim6"};
static const size_t dim_lens[VAR_RANK] = {7, 4, 2, 3, 5, 6, 9};

int
main(int argc, char **argv)
{
    /* mutually exclusive command line options */
	int option_group = 0;
	int option_deflate = 0;
	int option_unlimited = 0;
    /* file metadata */
    int mode = NC_CLOBBER;
    int ncid, grpid;
    int ivarid, fvarid;
    int ivar_dims[VAR_RANK];
    int fvar_dims[VAR_RANK];
    int ivar_data[NVALS];
    float fvar_data[NVALS];
    int r, i;
    char* file_name = FILE_NAME;
    int unlimid;

    /* Parse command line */
    if(argc >= 2) {
	file_name = argv[1];
    }
    if(argc >= 3) {
	if(strcmp(argv[2],"group")==0) {
	    option_group = 1;
	    mode |= NC_NETCDF4;
	} else if(strcmp(argv[2],"deflate")==0) {
	    option_deflate = 1;
	    mode |= NC_NETCDF4;
	} else if(strcmp(argv[2],"unlimited")==0) {
	    option_unlimited = 1;
	} else {
	    fprintf(stderr,"usage: tst_chunking [<filename> [group|deflate|unlimited]]\n");
	    exit(1);
	}
    }

    printf("*** Creating chunkable test file %s...\n", file_name);
    if(option_deflate)
	printf("\toption: deflate\n");
    else if(option_unlimited)
	printf("\toption: unlimited\n");
    else if(option_group)
	printf("\toption: group\n");

    if (nc_create(file_name, mode, &ncid)) LERR;
    for(r = 0; r < VAR_RANK; r++) {
	if (nc_def_dim(ncid, dim_names[r], dim_lens[r], &ivar_dims[r])) LERR;
	fvar_dims[VAR_RANK - 1 - r] = ivar_dims[r];
    }
    if(option_unlimited) {
	int udims[2];
	if (nc_def_dim(ncid, UNLIM_NAME, 0, &unlimid)) LERR;
	udims[0] = unlimid;
	udims[1] = ivar_dims[0];
	if (nc_def_var(ncid, IVAR_NAME, NC_INT, 2, udims, &ivarid)) LERR;
    } else {
	if (option_group) {
	    if (nc_def_grp(ncid, GRP_NAME, &grpid)) LERR;
	} else
	    grpid = ncid;
	if (nc_def_var(grpid, IVAR_NAME, NC_INT, VAR_RANK, ivar_dims, &ivarid)) LERR;
	if(option_deflate) {
	    if(nc_def_var_deflate(grpid,ivarid,NC_NOSHUFFLE, option_deflate, DEFLATE_LEVEL)) LERR;
	}
    }
    /* fvar is unchanged */
    if (nc_def_var(ncid, FVAR_NAME, NC_FLOAT, VAR_RANK, fvar_dims, &fvarid)) LERR;
    if (nc_enddef (ncid)) LERR;

    /* Fill in the data */
    if(option_unlimited) {
	int nvals = UNLIM_SIZE * dim_lens[0];
	size_t start[2] = {0,0};
	size_t count[2];
	for(i=0;i<nvals;i++) {
	    ivar_data[i] = i;	   
	}
	count[0] = UNLIM_SIZE;	
	count[1] = dim_lens[0];
	if (nc_put_vara(ncid, ivarid, start, count, ivar_data)) LERR;
    } else {
	for(i=0; i < NVALS; i++) {
	    ivar_data[i] = i;
	}
	if (nc_put_var(ncid, ivarid, ivar_data)) LERR;
    }
    /* fvar is unchanged */
    for(i=0; i < NVALS; i++) {
        fvar_data[i] = NVALS - i;
    }
    if (nc_put_var(ncid, fvarid, fvar_data)) LERR;

    if (nc_close(ncid)) LERR;

    SUMMARIZE_ERR;
    FINAL_RESULTS;
}
