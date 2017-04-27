/*
 * Demonstrate netcdf-4 rename bug.
 */

#include <netcdf.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* On error, prints line number and file of test program. */
#define ERR do {                                     \
fflush(stdout);                                      \
fprintf(stderr, "Unexpected result, %s, line: %d\n", \
	__FILE__, __LINE__);                         \
return 2;                                            \
} while (0)

#define FILE_NAME3 "tst_rnfix3.nc"
#define FILE_NAME4 "tst_rnfix4.nc"
#define ODIM_NAME "lat"		/* name for coord dim */
#define NDIM_NAME "tal"		/* new name for coord dim */
#define OVAR_NAME "lat"		/* name for coord var */
#define NVAR_NAME "tal"		/* new name for coord var */
#define OVAR2_NAME "rh"		/* name for non-coord var that uses coord dim */
#define VAR_RANK 1		/* all vars in this test are of same rank */
#define DIM_LEN 2		/* all dims in this test are of same len */

/* For renaming tests.  Create small test file of specified format
 * with a coordinate dimension, corresponding coordinate variable, and
 * a non-coordinate variable that uses the coordinate dimension.
 */
int
create_test_file(
    char *path,	/* filename */
    int format	/* NC_FORMAT_CLASSIC, NC_FORMAT_64BIT,
		   NC_FORMAT_NETCDF4, or NC_FORMAT_NETCDF4_CLASSIC */
    ) 
{
    int ncid, dimid, varid, var2id;
    int dims[VAR_RANK];
    int lats[DIM_LEN] = {-90, 90};
    float rh[DIM_LEN] = {0.25, 0.75};
    switch (format) {
    case (NC_FORMAT_CLASSIC):
	if (nc_create(path, 0, &ncid)) ERR;
	break;
    case (NC_FORMAT_64BIT_OFFSET):
	if (nc_create(path, NC_64BIT_OFFSET, &ncid)) ERR;
	break;
    case (NC_FORMAT_NETCDF4):
	if (nc_create(path, NC_NETCDF4, &ncid)) ERR;
	break;
    case(NC_FORMAT_NETCDF4_CLASSIC):
	if (nc_create(path, NC_NETCDF4 | NC_CLASSIC_MODEL, &ncid)) ERR;
	break;
    default:
	ERR;
	return NC_ENOTNC;
    }    
    if (nc_def_dim(ncid, ODIM_NAME, DIM_LEN, &dimid)) ERR;
    dims[0] = dimid;
    if (nc_def_var(ncid, OVAR_NAME, NC_INT, VAR_RANK, dims, &varid)) ERR;
    if (nc_def_var(ncid, OVAR2_NAME, NC_FLOAT, VAR_RANK, dims, &var2id)) ERR;
    if (nc_enddef(ncid)) ERR;	/* not necessary for netCDF-4 files */
    if (nc_put_var_int(ncid, varid, lats)) ERR;
    if (nc_put_var_float(ncid, var2id, rh)) ERR;
    if (nc_close(ncid)) ERR;
    return 0;
}
      
int
main(int argc, char **argv)
{
#define NUM_FORMATS 2
  int formats[NUM_FORMATS] = {NC_FORMAT_NETCDF4, NC_FORMAT_NETCDF4_CLASSIC};
  char *fmt_names[] = {"netCDF-4", "netCDF-4 classic model"};
  char *file_names[] = {FILE_NAME3, FILE_NAME4};
  int format;

  fprintf(stderr,"*** Testing netcdf rename bugs and fixes.\n");

  for(format = 0; format < NUM_FORMATS; format++)
    {
      int ncid, dimid, varid, var2id;
      int lats[DIM_LEN] = {-90, 90};
      int lats_in[DIM_LEN];
      float rh[DIM_LEN] = {0.25, 0.75};
      float rh_in[DIM_LEN];
      int ii;
      
      fprintf(stderr,"*** Test renaming coordinate variable and its dimension for %s...\n",
	     fmt_names[format]);
      if (create_test_file(file_names[format], formats[format])) ERR;
      if (nc_open(file_names[format], NC_WRITE, &ncid)) ERR;
      if (nc_inq_dimid(ncid, ODIM_NAME, &dimid)) ERR;
      if (nc_inq_varid(ncid, OVAR_NAME, &varid)) ERR;
      if (nc_inq_varid(ncid, OVAR2_NAME, &var2id)) ERR;
      if (nc_redef(ncid)) ERR; /* omitting this and nc_enddef call eliminates bug */
      if (nc_rename_dim(ncid, dimid, NDIM_NAME)) ERR;
      if (nc_rename_var(ncid, varid, NVAR_NAME)) ERR;
      if (nc_enddef(ncid)) ERR;
      if (nc_get_var_int(ncid, varid, lats_in)) ERR;
      for (ii = 0; ii < DIM_LEN; ii++) {
	if (lats_in[ii] != lats[ii])
	  fprintf(stderr, "\tlats_in[%d] is %d, should be %d\n", ii, lats_in[ii], lats[ii]);
      }
      if (nc_get_var_float(ncid, var2id, rh_in)) ERR;
      for (ii = 0; ii < DIM_LEN; ii++) {
	if (rh_in[ii] != rh[ii])
	  fprintf(stderr, "\trh_in[%d] is %g, should be %g\n", ii, rh_in[ii], rh[ii]);
      }
      if (nc_close(ncid)) ERR;

      fprintf(stderr,"*** Test renaming just coordinate variable for %s...\n",
	     fmt_names[format]);
      if (create_test_file(file_names[format], formats[format])) ERR;
      if (nc_open(file_names[format], NC_WRITE, &ncid)) ERR;
      if (nc_inq_dimid(ncid, ODIM_NAME, &dimid)) ERR;
      if (nc_inq_varid(ncid, OVAR_NAME, &varid)) ERR;
      if (nc_inq_varid(ncid, OVAR2_NAME, &var2id)) ERR;
      if (nc_redef(ncid)) ERR;  /* omitting this and nc_enddef call eliminates bug */
      /* if (nc_rename_dim(ncid, dimid, NDIM_NAME)) ERR; */
      if (nc_rename_var(ncid, varid, NVAR_NAME)) ERR;
      if (nc_enddef(ncid)) ERR;
      if (nc_get_var_int(ncid, varid, lats_in)) ERR;
      for (ii = 0; ii < DIM_LEN; ii++) {
	if (lats_in[ii] != lats[ii])
	  fprintf(stderr, "\tlats_in[%d] is %d, should be %d\n", ii, lats_in[ii], lats[ii]);
      }
      if (nc_get_var_float(ncid, var2id, rh_in)) ERR;
      for (ii = 0; ii < DIM_LEN; ii++) {
	if (rh_in[ii] != rh[ii])
	  fprintf(stderr, "\trh_in[%d] is %g, should be %g\n", ii, rh_in[ii], rh[ii]);
      }
      if (nc_close(ncid)) ERR;

      
      fprintf(stderr,"*** Test renaming just coordinate dimension for %s...\n",
	     fmt_names[format]);
      if (create_test_file(file_names[format], formats[format])) ERR;
      if (nc_open(file_names[format], NC_WRITE, &ncid)) ERR;
      if (nc_inq_dimid(ncid, ODIM_NAME, &dimid)) ERR;
      if (nc_inq_varid(ncid, OVAR_NAME, &varid)) ERR;
      if (nc_inq_varid(ncid, OVAR2_NAME, &var2id)) ERR;
      if (nc_redef(ncid)) ERR; /* omitting this and nc_enddef call eliminates bug */
      if (nc_rename_dim(ncid, dimid, NDIM_NAME)) ERR;
      /* if (nc_rename_var(ncid, varid, NVAR_NAME)) ERR; */
      if (nc_enddef(ncid)) ERR;
      if (nc_get_var_int(ncid, varid, lats_in)) ERR;
      for (ii = 0; ii < DIM_LEN; ii++) {
	if (lats_in[ii] != lats[ii])
	  fprintf(stderr, "\tlats_in[%d] is %d, should be %d\n", ii, lats_in[ii], lats[ii]);
      }
      if (nc_get_var_float(ncid, var2id, rh_in)) ERR;
      for (ii = 0; ii < DIM_LEN; ii++) {
	if (rh_in[ii] != rh[ii])
	  fprintf(stderr, "\trh_in[%d] is %g, should be %g\n", ii, rh_in[ii], rh[ii]);
      }
      if (nc_close(ncid)) ERR;

      if (formats[format] == NC_FORMAT_NETCDF4) {
          printf("*** Test renaming attribute in sub-group for %s...\n",
                fmt_names[format]);
          {
#define DIMNAME "lon"
#define VARNAME "lon"
#define G1_VARNAME "lon"
#define OLD_NAME "units"
#define NEW_NAME "new_units"
#define CONTENTS "degrees_east"
#define RANK_lon 1
#define GRP_NAME "g1"
#define RANK_g1_lon 1

            /* IDs of file, groups, dimensions, variables, attributes */
            int ncid, g1_grp, lon_dim, lon_var, g1_lon_var, units_att;
            size_t lon_len = 4;
            char *data_in;

            /* variable shapes */
            int lon_dims[RANK_lon];
            int g1_lon_dims[RANK_g1_lon];

            if (!(data_in = malloc(strlen(CONTENTS) + 1))) ERR;

            /* Create test file */
            if (nc_create(file_names[format], NC_NETCDF4 | NC_CLOBBER, &ncid)) ERR;
            /* Create subgroup and outer dimension */
            if (nc_def_grp(ncid, GRP_NAME, &g1_grp)) ERR;
            if (nc_def_dim(ncid, DIMNAME, lon_len, &lon_dim)) ERR;
            /* Create outer variable and subgroup variable */
            lon_dims[0] = lon_dim;
            if (nc_def_var(ncid, VARNAME, NC_FLOAT, RANK_lon, lon_dims, &lon_var)) ERR;
            g1_lon_dims[0] = lon_dim;
            if (nc_def_var(g1_grp, G1_VARNAME, NC_FLOAT, RANK_g1_lon, g1_lon_dims, &g1_lon_var)) ERR;
            /* assign per-variable attributes */
            if (nc_put_att_text(ncid, lon_var, OLD_NAME, strlen(CONTENTS), CONTENTS)) ERR;
            if (nc_put_att_text(g1_grp, g1_lon_var, OLD_NAME, strlen(CONTENTS), CONTENTS)) ERR;
            if (nc_enddef (ncid)) ERR;
            /* write variable data */
            {
              float lon_data[4] = {0, 90, 180, 270};
              size_t start[] = {0};
              size_t count[] = {4};
              if (nc_put_vara(ncid, lon_var, start, count, lon_data)) ERR;
            }
            {
              float g1_lon_data[4] = {0, 90, 180, 270};
              size_t start[] = {0};
              size_t count[] = {4};
              if (nc_put_vara(g1_grp, g1_lon_var, start, count, g1_lon_data)) ERR;
            }
            if (nc_close(ncid)) ERR;

            /* reopen the file and rename the attribute in the subgroup */
            if (nc_open(file_names[format], NC_WRITE, &ncid)) ERR;
            if (nc_inq_grp_ncid(ncid, GRP_NAME, &g1_grp)) ERR;
            if (nc_inq_varid(g1_grp,VARNAME,&g1_lon_var)) ERR;
            if (nc_rename_att(g1_grp, g1_lon_var, OLD_NAME, NEW_NAME)) ERR;
            if (nc_close(ncid)) ERR;

            /* reopen the file again and see if renamed attribute exists and
               has expected value */
            {
              nc_type att_type;
              size_t att_len;

              if (nc_open(file_names[format], NC_NOWRITE, &ncid)) ERR;
              if (nc_inq_grp_ncid(ncid, GRP_NAME, &g1_grp)) ERR;
              if (nc_inq_varid(g1_grp, VARNAME, &g1_lon_var)) ERR;
              if (nc_get_att_text(g1_grp, g1_lon_var, NEW_NAME, data_in)) ERR;
              if (strncmp(CONTENTS, data_in, strlen(CONTENTS))) ERR;
              if (nc_close(ncid)) ERR;
            }
            free(data_in);
          }
      }
    }

  return(0);
}
