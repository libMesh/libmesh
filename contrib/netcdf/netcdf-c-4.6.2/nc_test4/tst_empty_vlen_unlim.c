/* This is part of the netCDF package.  Copyright 2005 University
   Corporation for Atmospheric Research/Unidata See COPYRIGHT file for
   conditions of use.

   This program excersizes HDF5 variable length array code.

   $Id: tst_h_vl2.c,v 1.5 2010/06/01 15:34:52 ed Exp $
*/

/* This test was added to help diagnose the issue described at:
 *    * https://github.com/Unidata/netcdf-c/issues/221
 * This issue was originally reported by the python group at:
 *    * https://github.com/Unidata/netcdf4-python/issues/527
 */


#include <nc_tests.h>
#include "err_macros.h"
#include <hdf5.h>
#include <nc_logging.h>

#define FILE_NAME_UNLIM "tst_empty_vlen_unlim.nc"
#define FILE_NAME_LIM "tst_empty_vlen_lim.nc"
#define DIM_LEN_UNLIM NC_UNLIMITED
#define DIM_LEN_LIM 5
#define DIM_NAME "x"
#define VLEN_NAME "vltest"
#define VAR_NAME1 "v"
#define VAR_NAME2 "w"
#define ROW_COUNT 3
#define VLEN0 2
#define VLEN1 3
#define VLEN2 3

int main() {

  printf("Testing access to unset entries in VLEN variable, unlimited dimension\n");
  {
    int ncid, typeid, dimid, varid, varid2;
    nc_vlen_t data[ROW_COUNT];
    int stat;
    float *dat0, *dat1, *dat2;
    float *data2;
    size_t startp[3] = {0,0,0};
    size_t countp[3] = {VLEN0,VLEN1,VLEN2};
    size_t startp2[1] = {0};
    size_t countp2[1] = {VLEN2};
    int i = 0;
    /* Create File */
    printf("\t* Creating File:\tnc_create()\n");
    if (nc_create(FILE_NAME_UNLIM, NC_NETCDF4 | NC_CLOBBER, &ncid)) ERR;

    /* Set fill mode */
    //printf("\t* Setting fill mode:\tnc_set_fill()\n");
    //if(nc_set_fill(ncid,NC_FILL,NULL)) ERR;

    /* Create Dimension */
    printf("\t* Defining Unlimited Dimension:\tnc_def_dim()\n");
    if (nc_def_dim(ncid, DIM_NAME, DIM_LEN_UNLIM, &dimid)) ERR;

    /* Create ragged array type. */
    printf("\t* Creating Ragged Array type:\tnc_def_vlen().\n");
    if (nc_def_vlen(ncid, VLEN_NAME, NC_FLOAT, &typeid)) ERR;

    /* Create a variable of typeid. */
    printf("\t* Creating Variable using Ragged Arrayt Type:\tnc_def_var().\n");
    if (nc_def_var(ncid, VAR_NAME1, typeid, 1, &dimid, &varid)) ERR;

    /* Create a variable of type float. */
    printf("\t* Creating secondary Variable using NC_FLOAT:\tnc_def_var().\n");
    if (nc_def_var(ncid, VAR_NAME2, NC_FLOAT, 1, &dimid, &varid2)) ERR;

    /* End define mode. */
    printf("\t* Ending define mode:\tnc_enddef().\n");

    /* Write out data for w */
    printf("\t* Creating float data for secondary variable.\n");
    data2 = (float*)malloc(sizeof(float) * VLEN2);
    for(i = 0; i < VLEN2; i++) {
      data2[i] = (float)i;
    }

    printf("\t* Putting data in secondary variable:\tnc_put_vara().\n");
    if (nc_put_vara(ncid,varid2,startp2,countp2,data2)) ERR;
    free(data2);

    /***********/
    /* Actually unnecessary to recreate the issue. */
    /***********/


    /* Write out varying-length data for v[0] and v[1]. Leave v[2] empty. */

    dat0 = (float*)malloc(VLEN0 * sizeof(float));
    for(i = 0; i < VLEN0; i++) {
      dat0[i] = (float)i;
    }
    dat1 = (float*)malloc(VLEN1 * sizeof(float));
    for(i = 0; i < VLEN1; i++) {
      dat1[i] = (float)i;
    }
    dat2 = (float*)malloc(VLEN2 * sizeof(float));
    for(i = 0; i < VLEN2; i++) {
      dat2[i] = (float)i;
    }

    data[0].p = dat0;
    data[0].len = VLEN0;

    data[1].p = dat1;
    data[1].len = VLEN1;

    data[2].p = dat2;
    data[2].len = VLEN2;

    printf("\t* Putting data in VLEN variable:\tnc_put_vara().\n");
    stat = nc_put_vara(ncid,varid,startp,countp,data);
    if(stat) ERR;

    /* Close File. */
    printf("\t* Closing file:\tnc_close().\n");
    if ((stat = nc_close(ncid))) ERR;

    free(dat0);
    free(dat1);
    free(dat2);
  }

  printf("Testing access to unset entries in VLEN variable, unlimit dimension\n");
  {
    int ncid, typeid, dimid, varid, varid2;
    nc_vlen_t data[ROW_COUNT];
    int stat;
    float *dat0, *dat1, *dat2;
    float *data2;
    size_t startp[3] = {0,0,0};
    size_t countp[3] = {VLEN0,VLEN1,VLEN2};
    size_t startp2[1] = {0};
    size_t countp2[1] = {VLEN2};
    int i = 0;
    /* Create File */
    printf("\t* Creating File:\tnc_create()\n");
    if (nc_create(FILE_NAME_LIM, NC_NETCDF4 | NC_CLOBBER, &ncid)) ERR;

    /* Set fill mode */
    //printf("\t* Setting fill mode:\tnc_set_fill()\n");
    //if(nc_set_fill(ncid,NC_FILL,NULL)) ERR;

    /* Create Dimension */
    printf("\t* Defining Unlimited Dimension:\tnc_def_dim()\n");
    if (nc_def_dim(ncid, DIM_NAME, DIM_LEN_LIM, &dimid)) ERR;

    /* Create ragged array type. */
    printf("\t* Creating Ragged Array type:\tnc_def_vlen().\n");
    if (nc_def_vlen(ncid, VLEN_NAME, NC_FLOAT, &typeid)) ERR;

    /* Create a variable of typeid. */
    printf("\t* Creating Variable using Ragged Arrayt Type:\tnc_def_var().\n");
    if (nc_def_var(ncid, VAR_NAME1, typeid, 1, &dimid, &varid)) ERR;

    /* Create a variable of type float. */
    printf("\t* Creating secondary Variable using NC_FLOAT:\tnc_def_var().\n");
    if (nc_def_var(ncid, VAR_NAME2, NC_FLOAT, 1, &dimid, &varid2)) ERR;

    /* End define mode. */
    printf("\t* Ending define mode:\tnc_enddef().\n");

    /* Write out data for w */
    printf("\t* Creating float data for secondary variable.\n");
    data2 = (float*)malloc(sizeof(float) * VLEN2);
    for(i = 0; i < VLEN2; i++) {
      data2[i] = (float)i;
    }

    printf("\t* Putting data in secondary variable:\tnc_put_vara().\n");
    if (nc_put_vara(ncid,varid2,startp2,countp2,data2)) ERR;
    free(data2);

    /***********/
    /* Actually unnecessary to recreate the issue. */
    /***********/


    /* Write out varying-length data for v[0] and v[1]. Leave v[2] empty. */

    dat0 = (float*)malloc(VLEN0 * sizeof(float));
    for(i = 0; i < VLEN0; i++) {
      dat0[i] = (float)i;
    }
    dat1 = (float*)malloc(VLEN1 * sizeof(float));
    for(i = 0; i < VLEN1; i++) {
      dat1[i] = (float)i;
    }
    dat2 = (float*)malloc(VLEN2 * sizeof(float));
    for(i = 0; i < VLEN2; i++) {
      dat2[i] = (float)i;
    }

    data[0].p = dat0;
    data[0].len = VLEN0;

    data[1].p = dat1;
    data[1].len = VLEN1;

    data[2].p = dat2;
    data[2].len = VLEN2;

    printf("\t* Putting data in VLEN variable:\tnc_put_vara().\n");
    stat = nc_put_vara(ncid,varid,startp,countp,data);
    if(stat) ERR;

    /* Close File. */
    printf("\t* Closing file:\tnc_close().\n");
    if ((stat = nc_close(ncid))) ERR;
    free(dat0);
    free(dat1);
    free(dat2);
  }

  SUMMARIZE_ERR;
  FINAL_RESULTS;

}
