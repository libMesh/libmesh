/* Copyright 2010 University Corporation for Atmospheric Research. See
   COPYRIGHT file for copying and redistribution conditions.

   This test program is only built if netCDF-4 is disabled. It tests
   the netCDF-3 version of nc_inq_type().

 $Id: t_type.c,v 2.3 2010/05/26 21:43:33 dmh Exp $ 
 */

#include <config.h>
#include <nc_tests.h>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <netcdf.h>

int
main(int ac, char *av[])
{
   int ncid = 0;
/*   int t; */
   size_t size_in;
   char name_in[NC_MAX_NAME + 1];

   printf("\n*** Testing netCDF classic version of nc_inq_type...");
   if (nc_inq_type(ncid, 0, name_in, &size_in) != NC_EBADTYPE) ERR;
   if (nc_inq_type(ncid, NC_STRING + 1, name_in, &size_in) != NC_EBADTYPE) ERR;
   if (nc_inq_type(ncid, NC_BYTE, name_in, &size_in)) ERR;
   if (strcmp(name_in, "byte") || size_in != 1) ERR;
   if (nc_inq_type(ncid, NC_CHAR, name_in, &size_in)) ERR;
   if (strcmp(name_in, "char") || size_in != 1) ERR;
   if (nc_inq_type(ncid, NC_SHORT, name_in, &size_in)) ERR;
   if (strcmp(name_in, "short") || size_in != 2) ERR;
   if (nc_inq_type(ncid, NC_INT, name_in, &size_in)) ERR;
   if (strcmp(name_in, "int") || size_in != 4) ERR;
   if (nc_inq_type(ncid, NC_FLOAT, name_in, &size_in)) ERR;
   if (strcmp(name_in, "float") || size_in != 4) ERR;
   if (nc_inq_type(ncid, NC_DOUBLE, name_in, &size_in)) ERR;
   if (strcmp(name_in, "double") || size_in != 8) ERR;
   SUMMARIZE_ERR;
   FINAL_RESULTS;
}
