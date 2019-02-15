/* This is part of the netCDF package.
   Copyright 2006 University Corporation for Atmospheric Research/Unidata.
   See COPYRIGHT file for conditions of use.

   This test was provided by Jeff Whitaker as an example of a bug,
   specifically a segfault when re-writing an NC_CHAR attribute as
   an NC_STRING attribute.

   See https://github.com/Unidata/netcdf-c/issues/149

   $Id$
*/

#include <netcdf.h>
#include <config.h>
#include <nc_tests.h>
#include "err_macros.h"
#include <string.h>

#define FILE_NAME "tst_atts_string_rewrite.nc"

int main() {
   int dataset_id;
   const char *attstring[1] = {"bar"};
   int res = 0;
   printf("\n*** Testing overwriting text attribute with string attribute.\n");
   printf("\n***Creating file...");

   res = nc_create(FILE_NAME, NC_NETCDF4, &dataset_id); if(res) ERR;
   printf("Success\n");

   printf("Creating global attribute with nc_put_att_text...");
   res = nc_put_att_text(dataset_id, NC_GLOBAL, "foo", 3, "bar");
   printf("Success\n");

   printf("Overwriting global attribute with nc_put_att_string...");
   res = nc_put_att_string(dataset_id, NC_GLOBAL, "foo", 1, attstring);
   printf("Success\n");

   printf("Closing file...");
   res = nc_close(dataset_id);
   printf("Success\n");

   printf("Test Finished.\n");
   SUMMARIZE_ERR;
   FINAL_RESULTS;

}
