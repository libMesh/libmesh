/* This is part of the netCDF package.  Copyright 2005-2018,
   University Corporation for Atmospheric Research/Unidata. See
   COPYRIGHT file for conditions of use.

   Test that NetCDF-4 can read HDF4 files.
   Ed Hartnett, Ward Fisher
*/
#include <config.h>
#include <nc_tests.h>
#include "err_macros.h"
#include <mfhdf.h>
#include <netcdf_f.h>

#define FILE_NAME "tst_hdf4_extra.h4"

#define PRES_NAME "pres"
#define LAT_LEN 3
#define LON_LEN 2
#define NDIMS2 2
#define ATT_NAME "Caesar"
#define NAME_DUMB "Bozo"

int
create_hdf4_file()
{
   int32 sd_id, sds_id;
   int32 dim_size[NDIMS2] = {LAT_LEN, LON_LEN};
   int32 start[NDIMS2] = {0, 0}, edge[NDIMS2] = {LAT_LEN, LON_LEN};
   int data_out[LAT_LEN][LON_LEN];
   int test_val = 42;
   int i, j;
   int count = 0;

   /* Create some data. */
   for (i = 0; i < LAT_LEN; i++)
      for (j = 0; j < LON_LEN; j++)
         data_out[i][j] = count++;

   /* Create a file with one SDS, containing our phony data. */
   sd_id = SDstart(FILE_NAME, DFACC_CREATE);
   sds_id = SDcreate(sd_id, PRES_NAME, DFNT_INT32, NDIMS2, dim_size);
   if (SDwritedata(sds_id, start, NULL, edge, (void *)data_out)) ERR;

   /* Add a global attribute. */
   if (SDsetattr(sd_id, ATT_NAME, DFNT_INT32, 1, &test_val)) ERR;

   /* Shut down. */
   if (SDendaccess(sds_id)) ERR;
   if (SDend(sd_id)) ERR;

   return 0;
}

int
main(int argc, char **argv)
{
   printf("\n*** Testing HDF4/NetCDF-4 interoperability extra stuff...\n");

   /* Create our test file. */
   if (create_hdf4_file()) ERR;

   printf("*** testing data conversion...");
   {
      int ncid;
      size_t start[NDIMS2] = {0, 0}, count[NDIMS2] = {LAT_LEN, LON_LEN};
      int data_int[LAT_LEN * LON_LEN];
      short data_short[LAT_LEN * LON_LEN];
      int data_int2[LAT_LEN * LON_LEN];
      float data_float[LAT_LEN * LON_LEN];
      double data_double[LAT_LEN * LON_LEN];
      int i = 0;

      /* Open HDF4 file with netCDF. */
      if (nc_open(FILE_NAME, 0, &ncid)) ERR;

      /* These won't work. */
      if (nc_get_vara_int(ncid, 0, NULL, count, data_int) != NC_EINVALCOORDS) ERR;
      if (nc_get_vara_int(ncid + TEST_VAL_42, 0, start, count, data_int) != NC_EBADID) ERR;

      /* Read data as short. */
      if (nc_get_vara_short(ncid, 0, start, count, data_short)) ERR;
      for (i = 0; i < LAT_LEN * LON_LEN; i++)
         if (data_short[i] != (short)i) ERR;

      /* Read data as int. */
      if (nc_get_vara_int(ncid, 0, start, count, data_int)) ERR;
      for (i = 0; i < LAT_LEN * LON_LEN; i++)
         if (data_int[i] != i) ERR;

      /* NULL count is treated as meaing entire variable. */
      if (nc_get_vara_int(ncid, 0, start, NULL, data_int2)) ERR;
      for (i = 0; i < LAT_LEN * LON_LEN; i++)
         if (data_int2[i] != i) ERR;

      /* Read data as float. */
      if (nc_get_vara_float(ncid, 0, start, count, data_float)) ERR;
      for (i = 0; i < LAT_LEN * LON_LEN; i++)
         if (data_float[i] != (float)i) ERR;

      /* Read data as double. */
      if (nc_get_vara_double(ncid, 0, start, count, data_double)) ERR;
      for (i = 0; i < LAT_LEN * LON_LEN; i++)
         if (data_double[i] != (double)i) ERR;

      /* Close the file. */
      if (nc_close(ncid)) ERR;
   }
   SUMMARIZE_ERR;
   printf("*** testing bad parameters, read-only writes, and abort...");
   {
      int ncid;
      int ndims, nvars, ngatts, unlimdimid;
      size_t start[NDIMS2] = {0, 0}, count[NDIMS2] = {1, 1};
      int test_val;

      /* These will not work. */
      if (nc_open(FILE_NAME, NC_MMAP, &ncid) != NC_EINVAL) ERR;
      if (nc_open(FILE_NAME, NC_64BIT_OFFSET, &ncid) != NC_EINVAL) ERR;
      if (nc_open(FILE_NAME, NC_DISKLESS, &ncid) != NC_EINVAL) ERR;

      /* Now open with netCDF. */
      if (nc_open(FILE_NAME, 0, &ncid)) ERR;

      /* Check it out. */
      if (nc_inq(ncid, &ndims, &nvars, &ngatts, &unlimdimid)) ERR;
      if (ndims != 2 || nvars != 1 || ngatts != 1 || unlimdimid != -1) ERR;
      if (nc_inq(ncid, NULL, NULL, NULL, NULL)) ERR;
      if (nc_inq(ncid + TEST_VAL_42, NULL, NULL, NULL, NULL) != NC_EBADID) ERR;

      /* These only work for netCDF-3 files. */
      if (nc_set_base_pe(ncid, 0) != NC_ENOTNC3) ERR;
      if (nc_inq_base_pe(ncid, NULL) != NC_ENOTNC3) ERR;

      /* Attempt to write. */
      if (nc_rename_att(ncid, NC_GLOBAL, ATT_NAME, NAME_DUMB) != NC_EPERM) ERR;
      if (nc_del_att(ncid, NC_GLOBAL, ATT_NAME) != NC_EPERM) ERR;
      if (nc_put_att_int(ncid, NC_GLOBAL, NAME_DUMB, NC_INT, 0, NULL) != NC_EPERM) ERR;
      if (nc_def_dim(ncid, NAME_DUMB, 1, NULL) != NC_EPERM) ERR;
      if (nc_def_var(ncid, "hh", NC_INT, 0, NULL, NULL) != NC_EPERM) ERR;
      if (nc_def_var_fill(ncid, 0, 0, &test_val) != NC_EPERM) ERR;
      if (nc_rename_var(ncid, 0, NAME_DUMB) != NC_EPERM) ERR;
      if (nc_put_vara_int(ncid, 0, start, count, &test_val) != NC_EPERM) ERR;
      if (nc_set_fill(ncid, 0, NULL) != NC_EPERM) ERR;
      if (nc_rename_dim(ncid, 0, NULL) != NC_EPERM) ERR;

      /* These succeed but do nothing. */
      if (nc_enddef(ncid)) ERR;
      if (nc_sync(ncid)) ERR;

      /* These netcdf-4 operations are not supported. */
      if (nc_def_var_filter(ncid, 0, 0, 0, NULL) != NC_ENOTNC4) ERR;
      if (nc_def_var_fletcher32(ncid, 0, 0) != NC_ENOTNC4) ERR;
      if (nc_def_var_endian(ncid, 0, 0) != NC_ENOTNC4) ERR;
      if (nc_def_grp(ncid, NAME_DUMB, NULL) != NC_ENOTNC4) ERR;
      if (nc_rename_grp(ncid, NAME_DUMB) != NC_ENOTNC4) ERR;
      if (nc_def_compound(ncid, 1, NAME_DUMB, NULL) != NC_ENOTNC4) ERR;
      if (nc_insert_compound(ncid, 1, NAME_DUMB, 1, 1) != NC_ENOTNC4) ERR;
      if (nc_insert_array_compound(ncid, 1, NAME_DUMB, 1, 1, 1, NULL) != NC_ENOTNC4) ERR;
      if (nc_inq_compound_field(ncid, 1, 1, NULL, NULL, NULL, NULL, NULL) != NC_ENOTNC4) ERR;
      if (nc_inq_compound_fieldindex(ncid, 1, NULL, NULL) != NC_ENOTNC4) ERR;
      if (nc_def_opaque(ncid, 1, NULL, NULL) != NC_ENOTNC4) ERR;
      if (nc_def_vlen(ncid, NULL, 1, NULL) != NC_ENOTNC4) ERR;
      if (nc_def_enum(ncid, 1, NULL, NULL) != NC_ENOTNC4) ERR;
      if (nc_inq_enum_ident(ncid, 1, 1, NULL) != NC_ENOTNC4) ERR;
      if (nc_inq_enum_member(ncid, 1, 1, NULL, NULL) != NC_ENOTNC4) ERR;
      if (nc_insert_enum(ncid, 1, NULL, NULL) != NC_ENOTNC4) ERR;
      if (nc_put_vlen_element(ncid, 1, NULL, 1, NULL) != NC_ENOTNC4) ERR;
      if (nc_get_vlen_element(ncid, 1, NULL, NULL, NULL) != NC_ENOTNC4) ERR;
      if (nc_set_var_chunk_cache(ncid, 1, 1, 1, 1.0) != NC_ENOTNC4) ERR;
      if (nc_get_var_chunk_cache(ncid, 1, NULL, NULL, NULL) != NC_ENOTNC4) ERR;

      /* Abort is the same as nc_close, since HDF4 is read-only.) */
      if (nc_abort(ncid)) ERR;
   }
   SUMMARIZE_ERR;
   FINAL_RESULTS;
}
