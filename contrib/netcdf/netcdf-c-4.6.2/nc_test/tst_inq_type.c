/* This is part of the netCDF package. Copyright 2016 University
   Corporation for Atmospheric Research/Unidata See COPYRIGHT file for
   conditions of use. See www.unidata.ucar.edu for more info.

   Test nc_inq_type

   Added in support of https://github.com/Unidata/netcdf/issues/240

*/

#include <stdlib.h>
#include <string.h>

#include "config.h"
#include <nc_tests.h>
#include "err_macros.h"
#include <netcdf.h>

#ifdef USE_PNETCDF
#include <netcdf_par.h>
#endif

#define FILE_NAME "tst_inq_type.nc"

int test_type_should_fail(int ncid, int type, char* tstring) {

   printf("\t* Testing Type (Should Fail) %s:\t",tstring);
   if(!nc_inq_type(ncid,type,NULL,NULL)) ERR;
   else printf("expected failure.\n");

   return 0;
}

int test_type(int ncid, int type, char* tstring) {

   printf("\t* Testing Type %s:\t",tstring);
   if(nc_inq_type(ncid,type,NULL,NULL)) ERR;
   else printf("success.\n");

   return 0;
}

int main(int argc, char **argv) {

   int ncid=0;

   printf("\n* Testing nc_inq_type with netcdf-3\n");
   {
      if(nc_create(FILE_NAME,NC_CLOBBER,&ncid)) ERR;

      test_type(ncid, NC_BYTE,"NC_BYTE");
      test_type(ncid, NC_CHAR,"NC_CHAR");
      test_type(ncid, NC_SHORT,"NC_SHORT");
      test_type(ncid, NC_INT,"NC_INT");
      test_type(ncid, NC_LONG,"NC_LONG");
      test_type(ncid, NC_FLOAT,"NC_FLOAT");
      test_type(ncid, NC_DOUBLE,"NC_DOUBLE");

      /* Not Valid for Classic */
      /* Valid now, see https://github.com/Unidata/netcdf-c/issues/240 for more
	 information. The types are not valid for use in Classic,
	 but nc_inq_type should return valid info. */
      test_type(ncid, NC_UBYTE,"NC_UBYTE");
      test_type(ncid, NC_USHORT,"NC_USHORT");
      test_type(ncid, NC_UINT,"NC_UINT");
      test_type(ncid, NC_INT64,"NC_INT64");
      test_type(ncid, NC_UINT64,"NC_UINT64");
      test_type(ncid, NC_STRING,"NC_STRING");

      /* Invoke a true negative */
      test_type_should_fail(ncid, 9999, "NC_GARBAGE");
      test_type_should_fail(ncid, -1, "NC_GARBAGE_NEGATIVE");

      if(nc_close(ncid)) ERR;

      /* Reopen file to check that we can. */
      if (nc_create(FILE_NAME, 0, &ncid)) ERR;
      if (nc_close(ncid)) ERR;
      if (nc_create(FILE_NAME, NC_WRITE, &ncid)) ERR;
      if (nc_close(ncid)) ERR;
   }
   SUMMARIZE_ERR;

#ifdef ENABLE_CDF5
   printf("\n* Testing nc_inq_type with CDF5\n");
   {
      if(nc_create(FILE_NAME,NC_CLOBBER|NC_CDF5,&ncid)) ERR;

      test_type(ncid, NC_BYTE,"NC_BYTE");
      test_type(ncid, NC_CHAR,"NC_CHAR");
      test_type(ncid, NC_SHORT,"NC_SHORT");
      test_type(ncid, NC_INT,"NC_INT");
      test_type(ncid, NC_LONG,"NC_LONG");
      test_type(ncid, NC_FLOAT,"NC_FLOAT");
      test_type(ncid, NC_DOUBLE,"NC_DOUBLE");
      test_type(ncid, NC_UBYTE,"NC_UBYTE");
      test_type(ncid, NC_USHORT,"NC_USHORT");
      test_type(ncid, NC_UINT,"NC_UINT");
      test_type(ncid, NC_INT64,"NC_INT64");
      test_type(ncid, NC_UINT64,"NC_UINT64");
      test_type(ncid, NC_STRING,"NC_STRING");

      if(nc_close(ncid)) ERR;

      /* Reopen file to check that we can. */
      if (nc_open(FILE_NAME, NC_CDF5, &ncid)) ERR;
      if (nc_close(ncid)) ERR;
      if (nc_open(FILE_NAME, 0, &ncid)) ERR;
      if (nc_close(ncid)) ERR;
      if (nc_open(FILE_NAME, NC_WRITE, &ncid)) ERR;
      if (nc_close(ncid)) ERR;
   }
   SUMMARIZE_ERR;
#endif /* ENABLE_CDF5 */

#ifdef USE_NETCDF4
   printf("\n* Testing nc_inq_type with netcdf-4 + Classic Model\n");
   {
      if(nc_create(FILE_NAME,NC_CLOBBER|NC_NETCDF4|NC_CLASSIC_MODEL,&ncid)) ERR;

      test_type(ncid, NC_BYTE,"NC_BYTE");
      test_type(ncid, NC_CHAR,"NC_CHAR");
      test_type(ncid, NC_SHORT,"NC_SHORT");
      test_type(ncid, NC_INT,"NC_INT");
      test_type(ncid, NC_LONG,"NC_LONG");
      test_type(ncid, NC_FLOAT,"NC_FLOAT");
      test_type(ncid, NC_DOUBLE,"NC_DOUBLE");
      test_type(ncid, NC_UBYTE,"NC_UBYTE");
      test_type(ncid, NC_USHORT,"NC_USHORT");
      test_type(ncid, NC_UINT,"NC_UINT");
      test_type(ncid, NC_INT64,"NC_INT64");
      test_type(ncid, NC_UINT64,"NC_UINT64");
      test_type(ncid, NC_STRING,"NC_STRING");

      if(nc_close(ncid)) ERR;

      /* Re-open file to be sure we can. */
      if (nc_open(FILE_NAME, NC_NETCDF4|NC_CLASSIC_MODEL, &ncid)) ERR;
      if (nc_close(ncid)) ERR;
      if (nc_open(FILE_NAME, NC_NETCDF4, &ncid)) ERR;
      if (nc_close(ncid)) ERR;
      if (nc_open(FILE_NAME, NC_CLASSIC_MODEL, &ncid)) ERR;
      if (nc_close(ncid)) ERR;
      if (nc_open(FILE_NAME, 0, &ncid)) ERR;
      if (nc_close(ncid)) ERR;
      if (nc_open(FILE_NAME, NC_WRITE, &ncid)) ERR;
      if (nc_close(ncid)) ERR;
   }
   SUMMARIZE_ERR;

   printf("\n* Testing nc_inq_type with netcdf-4\n");
   {

      if(nc_create(FILE_NAME,NC_CLOBBER|NC_NETCDF4,&ncid)) ERR;

      test_type(ncid, NC_BYTE,"NC_BYTE");
      test_type(ncid, NC_CHAR,"NC_CHAR");
      test_type(ncid, NC_SHORT,"NC_SHORT");
      test_type(ncid, NC_INT,"NC_INT");
      test_type(ncid, NC_LONG,"NC_LONG");
      test_type(ncid, NC_FLOAT,"NC_FLOAT");
      test_type(ncid, NC_DOUBLE,"NC_DOUBLE");
      test_type(ncid, NC_UBYTE,"NC_UBYTE");
      test_type(ncid, NC_USHORT,"NC_USHORT");
      test_type(ncid, NC_UINT,"NC_UINT");
      test_type(ncid, NC_INT64,"NC_INT64");
      test_type(ncid, NC_UINT64,"NC_UINT64");
      test_type(ncid, NC_STRING,"NC_STRING");
      if(nc_close(ncid)) ERR;

      /* Re-open file to be sure we can. */
      if (nc_open(FILE_NAME, NC_NETCDF4, &ncid)) ERR;
      if (nc_close(ncid)) ERR;
      if (nc_open(FILE_NAME, 0, &ncid)) ERR;
      if (nc_close(ncid)) ERR;
      if (nc_open(FILE_NAME, NC_WRITE, &ncid)) ERR;
      if (nc_close(ncid)) ERR;
   }
   SUMMARIZE_ERR;

#endif // USE_NETCDF4

   printf("* Finished.\n");

   FINAL_RESULTS;
}

