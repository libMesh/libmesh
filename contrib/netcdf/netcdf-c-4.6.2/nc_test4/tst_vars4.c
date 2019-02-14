/* This is part of the netCDF package.
   Copyright 2005 University Corporation for Atmospheric Research/Unidata
   See COPYRIGHT file for conditions of use.

   Test netcdf-4 variables.
   Ed Hartnett
*/

#include <nc_tests.h>
#include "err_macros.h"

#define FILE_NAME "tst_vars4.nc"
#define NDIMS2 2
#define NUM_VARS 1
#define Y_NAME "y"
#define X_NAME "x"
#define VAR_NAME Y_NAME
#define XDIM_LEN 2
#define YDIM_LEN 5

int
main(int argc, char **argv)
{
   printf("\n*** Testing netcdf-4 variable functions, even more.\n");
   printf("**** testing Jeff's dimension problem...");
   {
      int varid, ncid, dims[NDIMS2], dims_in[NDIMS2];
      int ndims, nvars, ngatts, unlimdimid, natts;
      char name_in[NC_MAX_NAME + 1];
      nc_type type_in;
      size_t len_in;

      if (nc_create(FILE_NAME, NC_NETCDF4 | NC_CLOBBER, &ncid)) ERR;
      if (nc_def_dim(ncid, X_NAME, XDIM_LEN, &dims[0])) ERR;
      if (nc_def_dim(ncid, Y_NAME, YDIM_LEN, &dims[1])) ERR;
      if (nc_def_var(ncid, VAR_NAME, NC_FLOAT, 2, dims, &varid)) ERR;
      if (nc_inq(ncid, &ndims, &nvars, &ngatts, &unlimdimid)) ERR;
      if (nvars != NUM_VARS || ndims != NDIMS2 || ngatts != 0 || unlimdimid != -1) ERR;
      if (nc_inq_var(ncid, 0, name_in, &type_in, &ndims, dims_in, &natts)) ERR;
      if (strcmp(name_in, VAR_NAME) || type_in != NC_FLOAT || ndims != NDIMS2 ||
	  dims_in[0] != dims[0] || dims_in[1] != dims[1] || natts != 0) ERR;
      if (nc_inq_dim(ncid, 0, name_in, &len_in)) ERR;
      if (strcmp(name_in, X_NAME) || len_in != XDIM_LEN) ERR;
      if (nc_inq_dim(ncid, 1, name_in, &len_in)) ERR;
      if (strcmp(name_in, Y_NAME)) ERR;
      if (len_in != YDIM_LEN) ERR;
      if (nc_close(ncid)) ERR;

      /* Open the file and check. */
      if (nc_open(FILE_NAME, NC_WRITE, &ncid)) ERR;
      if (nc_inq(ncid, &ndims, &nvars, &ngatts, &unlimdimid)) ERR;
      if (nvars != NUM_VARS || ndims != NDIMS2 || ngatts != 0 || unlimdimid != -1) ERR;
      if (nc_inq_var(ncid, 0, name_in, &type_in, &ndims, dims_in, &natts)) ERR;
      if (strcmp(name_in, VAR_NAME) || type_in != NC_FLOAT || ndims != NDIMS2 ||
	  dims_in[0] != dims[0] || dims_in[1] != dims[1] || natts != 0) ERR;
      if (nc_inq_dim(ncid, 0, name_in, &len_in)) ERR;
      if (strcmp(name_in, X_NAME) || len_in != XDIM_LEN) ERR;
      if (nc_inq_dim(ncid, 1, name_in, &len_in)) ERR;
      if (strcmp(name_in, Y_NAME)) ERR;
      if (len_in != YDIM_LEN) ERR;
      if (nc_close(ncid)) ERR;
   }
   SUMMARIZE_ERR;
   printf("**** testing chunking turned on by fletcher...");
   {
      int varid, ncid, dims[NDIMS2];
      int storage_in;
      size_t chunksizes_in[NDIMS2];

      if (nc_create(FILE_NAME, NC_NETCDF4 | NC_CLOBBER, &ncid)) ERR;
      if (nc_def_dim(ncid, X_NAME, XDIM_LEN, &dims[0])) ERR;
      if (nc_def_dim(ncid, Y_NAME, YDIM_LEN, &dims[1])) ERR;
      if (nc_def_var(ncid, VAR_NAME, NC_FLOAT, 2, dims, &varid)) ERR;
      if (nc_def_var_fletcher32(ncid, varid, NC_FLETCHER32)) ERR;
      if (nc_inq_var_chunking(ncid, varid, &storage_in, chunksizes_in)) ERR;
      if (chunksizes_in[0] != XDIM_LEN || chunksizes_in[1] != YDIM_LEN) ERR;
      if (nc_close(ncid)) ERR;

      /* Open the file and check. */
      if (nc_open(FILE_NAME, NC_WRITE, &ncid)) ERR;
      if (nc_inq_var_chunking(ncid, varid, &storage_in, chunksizes_in)) ERR;
      if (chunksizes_in[0] != XDIM_LEN || chunksizes_in[1] != YDIM_LEN) ERR;
      if (nc_close(ncid)) ERR;
   }
   SUMMARIZE_ERR;
   printf("**** testing chunking turned on by shuffle...");
   {
      int varid, ncid, dims[NDIMS2];
      int storage_in;
      size_t chunksizes_in[NDIMS2];

      if (nc_create(FILE_NAME, NC_NETCDF4 | NC_CLOBBER, &ncid)) ERR;
      if (nc_def_dim(ncid, X_NAME, XDIM_LEN, &dims[0])) ERR;
      if (nc_def_dim(ncid, Y_NAME, YDIM_LEN, &dims[1])) ERR;
      if (nc_def_var(ncid, VAR_NAME, NC_FLOAT, 2, dims, &varid)) ERR;
      if (nc_def_var_deflate(ncid, varid, NC_SHUFFLE, 0, 0)) ERR;
      if (nc_inq_var_chunking(ncid, varid, &storage_in, chunksizes_in)) ERR;
      if (chunksizes_in[0] != XDIM_LEN || chunksizes_in[1] != YDIM_LEN) ERR;
      if (nc_close(ncid)) ERR;

      /* Open the file and check. */
      if (nc_open(FILE_NAME, NC_WRITE, &ncid)) ERR;
      if (nc_inq_var_chunking(ncid, varid, &storage_in, chunksizes_in)) ERR;
      if (chunksizes_in[0] != XDIM_LEN || chunksizes_in[1] != YDIM_LEN) ERR;
      if (nc_close(ncid)) ERR;
   }
   SUMMARIZE_ERR;
#define DIM_NAME "Distance_from_Mayo"
#define VAR_NAME_2 "Rocky_Road_to_Dublin"
#define NDIMS1 1
#define NUM_RECORDS 3
   printf("**** testing extending var along unlimited dim with no coord var...");
   {
      int varid, ncid, dimid;
      int ndims, nvars, natts, unlimdimid;
      size_t dim_len_in, index;
      int data = TEST_VAL_42;

      /* Create the test file with one var, one unlimited dim. */
      if (nc_create(FILE_NAME, NC_NETCDF4 | NC_CLOBBER, &ncid)) ERR;
      if (nc_def_dim(ncid, DIM_NAME, NC_UNLIMITED, &dimid)) ERR;
      if (nc_def_var(ncid, VAR_NAME_2, NC_INT, NDIMS1, &dimid, &varid)) ERR;

      /* Write some records. */
      for (index = 0; index < NUM_RECORDS; index++)
         if (nc_put_var1_int(ncid, varid, &index, &data)) ERR;
      if (nc_close(ncid)) ERR;

      /* Open the file and check. */
      if (nc_open(FILE_NAME, NC_WRITE, &ncid)) ERR;
      if (nc_inq(ncid, &ndims, &nvars, &natts, &unlimdimid)) ERR;
      if (ndims != 1 || nvars != 1 || natts != 0 || unlimdimid != 0) ERR;
      if (nc_inq_dim(ncid, dimid, NULL, &dim_len_in)) ERR;
      if (dim_len_in != NUM_RECORDS) ERR;

      /* Now add more records. */
      for (index = 3; index < NUM_RECORDS * 2; index++)
         if (nc_put_var1_int(ncid, varid, &index, &data)) ERR;
      if (nc_inq_dim(ncid, dimid, NULL, &dim_len_in)) ERR;

      if (dim_len_in != NUM_RECORDS * 2) ERR;

      /* Close the file. */
      if (nc_close(ncid)) ERR;
   }
   SUMMARIZE_ERR;
   FINAL_RESULTS;
}
