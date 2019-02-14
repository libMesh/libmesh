/* This is part of the netCDF package.
   Copyright 2016 University Corporation for Atmospheric Research/Unidata
   See COPYRIGHT file for conditions of use.

   Based on tst_fillbug.c

   Added in support of https://github.com/Unidata/netcdf-c/issues/239

   The issue in a nutshell is that if _FillValue is added to a file
   defined without it, previously existing attributes with
   preexisting variables are wiped out, only the _FillValue attribute
   exists.

*/

#include <nc_tests.h>
#include "err_macros.h"
#include <stdio.h>
#include <stdlib.h>
#include <netcdf.h>

#define FILENAME "tst_fill_attr_vanish.nc"
#define RANK_Time 1
#define RANK_P 3
#define LEN 4

#define ATTNAME "TextAttribute"
#define ATTVAL  "This is a text attribute used for testing."


/*! Main function for tst_fill_attr_vanish.c
 *
 */
int main()
{
  int ncid, dimids[RANK_P], time_id, p_id, test_id, status;
  int test_data[1] = {1};
  size_t test_start[1] = {0}, test_count[1] = {1};
  int test_fill_val[] = {5};
  double data[1] = {3.14159};
  size_t start[1] = {0}, count[1] = {1};

  printf("\n*** Testing for a netCDF-4 fill-value bug.\n");
  printf("*** Creating a file with no _FillValue defined. ***\n");

  /* Create a 3D test file. */
  if (nc_create(FILENAME, NC_CLOBBER|NC_NETCDF4, &ncid)) ERR;
  if (nc_set_fill(ncid, NC_NOFILL, NULL)) ERR;

  /* define dimensions */
  if (nc_def_dim(ncid, "Time", NC_UNLIMITED, &dimids[0])) ERR;
  if (nc_def_dim(ncid, "X", 4, &dimids[2])) ERR;
  if (nc_def_dim(ncid, "Y", 3, &dimids[1])) ERR;

  /* define variables */
  if (nc_def_var(ncid, "Time", NC_DOUBLE, 1, dimids, &time_id)) ERR;
  if (nc_def_var(ncid, "P", NC_FLOAT, RANK_P, dimids, &p_id)) ERR;
  if (nc_def_var(ncid, "Test", NC_INT, 1, &dimids[1], &test_id)) ERR;

  /* Add a _FillValue attribute */

  if (nc_put_att_text(ncid, test_id, ATTNAME, strlen(ATTVAL), ATTVAL)) ERR;
  /* Add a value to the test variable */
  if (nc_put_vara(ncid, test_id, test_start, test_count, test_data)) ERR;

  /* Add one record in coordinate variable. */
  if (nc_put_vara(ncid, time_id, start, count, data)) ERR;

  /* That's it! */
  if (nc_close(ncid)) ERR;

  /********************************************/

  /* Reopen the file, add a fillvalue attribute. */
  if (nc_open(FILENAME, NC_NOCLOBBER|NC_WRITE, &ncid)) ERR;
  if (nc_redef(ncid)) ERR;
  if (nc_inq_varid(ncid, "Test", &test_id)) ERR;

  /* Query existing attribute. */
  {
    char *attval = malloc(sizeof(char) * strlen(ATTVAL));
    printf("**** Checking that attribute still exists:\t");
    if(nc_get_att_text(ncid,test_id,ATTNAME,attval)) {printf("Fail\n"); ERR;}
    else {printf("%s\n",attval);}
    free(attval);

  }

  printf("**** Expecting NC_ELATEFILL when adding _FillValue attribute if variable exists.\n");
  status = nc_put_att_int(ncid, test_id, "_FillValue", NC_INT, 1, test_fill_val);
  if (status != NC_ELATEFILL) {
      fflush(stdout); /* Make sure our stdout is synced with stderr. */
      err++;
      fprintf(stderr, "Sorry! Expecting NC_ELATEFILL but got %s, at file %s line: %d\n",
              nc_strerror(status), __FILE__, __LINE__);
      return 2;
  }

  /* Query existing attribute. */
  {
    char *attval = malloc(sizeof(char) * strlen(ATTVAL));
    printf("**** Checking that attribute still exists, pre-write:\t");
    if(nc_get_att_text(ncid,test_id,ATTNAME,attval)) {printf("Fail\n"); ERR;}
    else {printf("%s\n",attval);}
    free(attval);

  }

  /* Close file again. */
  printf( "**** Saving, closing file.\n");
  if (nc_close(ncid)) ERR;
  /********************************************/
  printf( "*** Reopening file.\n");
  /* Reopen the file, checking that all attributes are preserved. */
  if (nc_open(FILENAME, NC_NOCLOBBER|NC_WRITE, &ncid)) ERR;
  if (nc_redef(ncid)) ERR;
  if (nc_inq_varid(ncid, "Test", &test_id)) ERR;

  /* Query existing attribute. */
  {
    char *attval = malloc(sizeof(char) * strlen(ATTVAL));
    printf("**** Checking that attribute still exists:\t");
    if(nc_get_att_text(ncid,test_id,ATTNAME,attval)) {printf("Fail\n"); ERR;}
    else {printf("%s\n",attval);}
    free(attval);

  }

  if (nc_close(ncid)) ERR;
  /********************************************/

  SUMMARIZE_ERR;

  FINAL_RESULTS;
  return 0;
}
