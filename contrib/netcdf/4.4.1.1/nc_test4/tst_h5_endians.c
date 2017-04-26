/*! Test for NCF-331. Added May 11, 2015.
 *  See the following links for more information:
 *
 *  o Issue on GitHub: https://github.com/Unidata/netcdf-c/issues/112
 *  o Issue in JIRA:   https://bugtracking.unidata.ucar.edu/browse/NCF-331
 *
 * Test contributed by Jeff Whitaker
 */

#include <string.h>
#include <netcdf.h>
#include <stdio.h>
#include <nc_tests.h>
#include "nc_logging.h"

#define FILE_NAME_NC "tst_h5_endians.nc"

#define NDIM 10
#define NLON 20
#define DIM_NAME "x"
#define DIM_LEN 4
#define GRP_NAME "grp"
#define LE_FLOAT_VARNAME "fl_le"
#define BE_FLOAT_VARNAME "fl_be"
#define LE_INT_VARNAME "int_le"
#define BE_INT_VARNAME "int_be"
#define LE_DBL_VARNAME "dbl_le"
#define BE_DBL_VARNAME "dbl_be"

int main() {

  int ncid, dimid;
  int le_float_varid;
  int be_float_varid;
  int le_int_varid;
  int be_int_varid;
  int le_dbl_varid;
  int be_dbl_varid;
  int ed;
  int failures = 0;
  int retval = 0;

  printf("* Checking that endianness is properly read from file.\n");
  printf("** Generating test files.\n");
  /*
   * 1. Create a netcdf file with endianness as desired.
   */
  {

    printf("*** Creating a file via netcdf API: %s.\n",FILE_NAME_NC);
    retval = nc_create(FILE_NAME_NC, NC_NETCDF4 | NC_CLOBBER, &ncid);

    retval = nc_def_dim(ncid, DIM_NAME, NDIM, &dimid);

    /* Little-Endian Float */
    retval = nc_def_var(ncid, LE_FLOAT_VARNAME, NC_FLOAT, 1, &dimid, &le_float_varid);
    retval = nc_def_var_endian(ncid, le_float_varid, NC_ENDIAN_LITTLE);

    /* Big-Endian Float */
    retval = nc_def_var(ncid, BE_FLOAT_VARNAME, NC_FLOAT, 1, &dimid, &be_float_varid);
    retval = nc_def_var_endian(ncid, be_float_varid, NC_ENDIAN_BIG);

    /* Little-Endian Int */
    retval = nc_def_var(ncid, LE_INT_VARNAME, NC_INT, 1, &dimid, &le_int_varid);
    retval = nc_def_var_endian(ncid, le_int_varid, NC_ENDIAN_LITTLE);

    /* Big-Endian Int */
    retval = nc_def_var(ncid, BE_INT_VARNAME, NC_INT, 1, &dimid, &be_int_varid);
    retval = nc_def_var_endian(ncid, be_int_varid, NC_ENDIAN_BIG);

    /* Little-Endian Double */
    retval = nc_def_var(ncid, LE_DBL_VARNAME, NC_DOUBLE, 1, &dimid, &le_dbl_varid);
    retval = nc_def_var_endian(ncid, le_dbl_varid, NC_ENDIAN_LITTLE);

    /* Big-Endian Double */
    retval = nc_def_var(ncid, BE_DBL_VARNAME, NC_DOUBLE, 1, &dimid, &be_dbl_varid);
    retval = nc_def_var_endian(ncid, be_dbl_varid, NC_ENDIAN_BIG);


    retval = nc_close(ncid);
  }

  /*
   * 2. Reopen netcdf-generated file, check to see if the endianness attribute
   *    exists.
   */
  printf("** Checking test files.\n");
  {
    ncid = 0;
    le_float_varid = 0;
    be_float_varid = 0;
    le_int_varid = 0;
    be_int_varid = 0;
    le_dbl_varid = 0;
    be_dbl_varid = 0;

    printf("*** %s\n",FILE_NAME_NC);
    retval = nc_open(FILE_NAME_NC, NC_NETCDF4 | NC_NOWRITE, &ncid);

    retval = nc_inq_varid(ncid,LE_FLOAT_VARNAME,&le_float_varid);
    retval = nc_inq_varid(ncid,BE_FLOAT_VARNAME,&be_float_varid);
    retval = nc_inq_varid(ncid,LE_INT_VARNAME,&le_int_varid);
    retval = nc_inq_varid(ncid,BE_INT_VARNAME,&be_int_varid);
    retval = nc_inq_varid(ncid,LE_DBL_VARNAME,&le_dbl_varid);
    retval = nc_inq_varid(ncid,BE_DBL_VARNAME,&be_dbl_varid);

    printf("\tLittle-Endian Float...\t");
    retval = nc_inq_var_endian(ncid,le_float_varid,&ed);
    if(ed == NC_ENDIAN_LITTLE) printf("passed\n"); else {printf("failed\n"); failures++;}

    printf("\tBig-Endian Float...\t");
    retval = nc_inq_var_endian(ncid,be_float_varid,&ed);
    if(ed == NC_ENDIAN_BIG) printf("passed\n"); else {printf("failed\n"); failures++;}

    printf("\tLittle-Endian Int...\t");
    retval = nc_inq_var_endian(ncid,le_int_varid,&ed);
    if(ed == NC_ENDIAN_LITTLE) printf("passed\n"); else {printf("failed\n"); failures++;}

    printf("\tBig-Endian Int...\t");
    retval = nc_inq_var_endian(ncid,be_int_varid,&ed);
    if(ed == NC_ENDIAN_BIG) printf("passed\n"); else {printf("failed\n"); failures++;}

    printf("\tLittle-Endian Double...\t");
    retval = nc_inq_var_endian(ncid,le_dbl_varid,&ed);
    if(ed == NC_ENDIAN_LITTLE) printf("passed\n"); else {printf("failed\n"); failures++;}

    printf("\tBig-Endian Double...\t");
    retval = nc_inq_var_endian(ncid,be_dbl_varid,&ed);
    if(ed == NC_ENDIAN_BIG) printf("passed\n"); else {printf("failed\n"); failures++;}

    retval = nc_close(ncid);
  }

  printf("** Failures Returned: [%d]\n",failures);
  return failures;
}
