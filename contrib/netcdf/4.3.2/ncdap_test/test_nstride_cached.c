/*
Report from Ansely Manke:
I've attached a file with a short c program that
demonstrates what I'm seeing. I don't know that we'd call
this a constraint, so that was a misleading description.

The thing I'm mimicking is a sequence where I read the 2-D
coordinate variables, make some decisions about the size of
a subset region, and then define index ranges and strides
for reading other variables.  Reading the full variable
happens correctly, but the read with strides winds up
messing up the data.  Here it shows up returning some valid
data and some zero's, but in other sequences of events it
seems to be actually mis-ordered.

I notice that the results are correct if the initial read of
the variable uses counts that are size-1. So that might be a
clue.  And the read is correct if I either close and reopen
the dataset after reading the original full variable, or if
I don't do that step of first reading the whole grid. It's
also fine with NetCDF 4.2.1.1.

-Ansley

Problem was two-fold:
1. struct Getvara has the dsttype field as OC_Type rather than nc_type
2. the dap odometer code dapnew_segment was incorrectly
   handling strides > 1: specifically, the stop position was incorrect.

*/

/* vars_whoi_test */
/* acm 4/2013 */
/* ansley.b.manke@noaa.gov */

/* test nc_get_vars_float with calls similar to Ferret calls */


/*linked with:

cc vars_whoi_test.c -g -o vars_whoi_test_4211 /usr/local/netcdf_4211/lib/libnetcdf.a /usr/local/hdf5_189/lib/libhdf5_hl.a /usr/local/hdf5_189/lib/libhdf5.a /usr/local/lib/libz.a -L/usr/lib64 -lc -lm -lcurl
cc vars_whoi_test.c -g -o vars_whoi_test /home/users/ansley/local/x86_nc43/lib/libnetcdf.a /usr/local/hdf5_189/lib/libhdf5_hl.a /usr/local/hdf5_189/lib/libhdf5.a /usr/local/lib/libz.a -L/usr/lib64 -lc -lm -lcurl

*/

/* Closing and reopening the dataset between the two reads fixes the
incorrect data return */

/* Setting the count to one less in the full data read also fixes the
incorrect data return */

#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include "netcdf.h"

static int verbose = 0;

static char* URL="http://geoport.whoi.edu/thredds/dodsC/coawst_4/use/fmrc/coawst_4_use_best.ncd";

int
main()
{

    int ncid;
    int varid;
    int i;
    int ncstatus;
    size_t start[5], count[5];
    ptrdiff_t stride[5], tmp_ptrdiff_t;
    int pass = 1;

    int idim, ndim;
    //float dat[301060];
    float *dat = (float*)malloc(sizeof(float)*301060);
    float sdat[10];

    for (idim=0; idim<5; idim++) {
        start[idim] = 0;
        count[idim] = 1;
        stride[idim] = 1;
    }

    ndim=2;

    printf(" \n");
    printf("********************\n");
    printf("open URL %s\n",URL);
    printf(" \n");

    ncstatus = nc_open(URL, NC_NOWRITE, &ncid);

    if(ncstatus != NC_NOERR) {
	fprintf(stderr,"Could not open: %s; server may be down; test ignored\n",URL);
	exit(0);
    }

    ncstatus = nc_inq_varid(ncid, "lon_rho", &varid);

    ndim=2;


if(verbose) {
    printf(" \n");
    printf("********************\n");
    printf("Read lon_rho data w/o strides\n");
    printf(" \n");
}

    start[0] = 0;
    start[1] = 0;
    start[2] = 0;
    start[3] = 0;

    count[0] =  336;
    count[1] =  896;

    stride[0] = 1;
    stride[1] = 1;

if(verbose) {
    for (idim=0; idim<ndim; idim++)
	printf("start[%1d]=%3lu count[%1d]=%3lu stride[%1d]=%3lu\n",
		idim,start[idim],idim,count[idim],idim,stride[idim]);
}

    ncstatus = nc_get_vars_float (ncid, varid, start, count, stride, (float*) dat);

if(verbose) {
    printf(" \n");
    printf("********************\n");
    printf("Print some of lon_rho\n");
    printf(" \n");

    for (i=0; i<10; i++)
        printf("lon_rho[%d] = %f\n",i,dat[i]);
    printf(" \n");

    for (i=301045; i<301055; i++)
        printf("lon_rho[%d] = %f\n",i,dat[i]);
    printf(" \n");
}

    memset((void*)dat,0,sizeof(dat));

    /* Read a second variable */

    ncstatus = nc_inq_varid(ncid, "lat_rho", &varid);

    ndim=2;

if(verbose) {
    printf(" \n");
    printf("********************\n");
    printf("Read lat_rho data w/o strides\n");
    printf(" \n");
}

    start[0] = 0;
    start[1] = 0;
    start[2] = 0;
    start[3] = 0;

    count[0] =  336;
    count[1] =  896;

    stride[0] = 1;
    stride[1] = 1;

if(verbose) {
    for (idim=0; idim<ndim; idim++)
        printf("start[%d]=%3lu count[%d]=%3lu stride[%d]=%3lu\n",
		idim, start[idim], idim, count[idim], idim, stride[idim]);
}

    ncstatus = nc_get_vars_float (ncid, varid, start, count, stride, (float*) dat);

if(verbose) {
    printf(" \n");
    printf("********************\n");
    printf("Print some of lat_rho\n");
    printf(" \n");
    for (i=0; i<10; i++)
        printf("lat_rho[%d] = %f\n",i,dat[i]);
    printf(" \n");

    printf(" \n");
    for (i=301045; i<301055; i++)
        printf("lon_rho[%d] = %f\n",i,dat[i]);
    printf(" \n");
}

    memset((void*)dat,0,sizeof(dat));

    /* close and reopen the dataset, then the below read is correct */

if(verbose) {
    printf(" \n");
    printf("********************\n");
    printf("Close and reopen the dataset\n");
}

    ncstatus = nc_close (ncid);
    ncstatus = nc_open(URL, NC_NOWRITE, &ncid);

    /*  ----------------------------------------------------- */
    /* Read a subset of the data with strides */

    ncstatus = nc_inq_varid(ncid, "lon_rho", &varid);

if(verbose) {
    printf(" \n");
    printf("********************\n");
    printf("Read a subset of lon_rho data with strides\n");
    printf(" \n");
}

    start[0] = 250;
    start[1] = 704;

    count[0] =  5;
    count[1] =  2;

    stride[0] = 2;
    stride[1] = 4;

if(verbose) {
    for (idim=0; idim<ndim; idim++)
	printf("start[%1d]=%3lu count[%1d]=%3lu stride[%1d]=%3lu\n",
		idim,start[idim],idim,count[idim],idim,stride[idim]);
}

    memset((void*)sdat,0,sizeof(sdat));
    ncstatus = nc_get_vars_float (ncid, varid, start, count, stride,  (float*) sdat);

    printf("status = %d\n", ncstatus);

    /* Verify that all read values are 67 <= n < 68 */
    for (i=0; i<10; i++) {
	if(!(sdat[i] <= -67 && sdat[i] > -68)) {
	    printf("lon_rho[%d] = %f\n",i,sdat[i]);
	    pass = 0;
	}
    }
if(verbose) {
    printf(" \n");
    printf("********************\n");
    printf("Print  values read. They should all be -67.xxxx \n");
    printf(" \n");

    for (i=0; i<10; i++)
        printf("lon_rho[%d] = %f\n",i,sdat[i]);
}

    ncstatus = nc_close (ncid);

    if(!pass) {
	printf("*** FAIL: lon_rho value out of range.\n");
	exit(1);
    }

    printf("*** PASS\n");
    free(dat);
    exit(0);

}
