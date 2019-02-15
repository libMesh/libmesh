/* This is part of the netCDF package. Copyright 2005-2011, University
   Corporation for Atmospheric Research/Unidata. See COPYRIGHT file
   for conditions of use.

   Test that NetCDF-4 can read a bunch of HDF4 files pulled in from
   the FTP site.

   @author Ed Hartnett
*/

#include <config.h>
#include <nc_tests.h>
#include "err_macros.h"
#include <mfhdf.h>

#define FILE_NAME "tst_interops3.h4"

/* Function to test nc_inq_format(). */
int
check_inq_format(int ncid, int expected_format, int expected_extended_format, int expected_mode)
{
   int format;
   int extended_format;
   int mode;
   
   if (nc_inq_format(ncid + 66000, NULL) != NC_EBADID) ERR;
   if (nc_inq_format(ncid, NULL)) ERR;
   if (nc_inq_format(ncid, &format)) ERR;
   if (format != expected_format) {
      printf("format %d expected_format %d\n", format, expected_format);      
      ERR;
   }
   if (nc_inq_format_extended(ncid + 66000, &extended_format, &mode) != NC_EBADID) ERR;
   {
      int mode;
      if (nc_inq_format_extended(ncid, NULL, &mode)) ERR;
      if (mode != expected_mode) {
         printf("expected_mode %x mode %x\n", expected_mode, mode);
         //ERR;
      }
   }
   {
      int extended_format;
      if (nc_inq_format_extended(ncid, &extended_format, NULL)) ERR;
      if (extended_format != expected_extended_format) ERR;
   }

   if (nc_inq_format_extended(ncid, &extended_format, &mode)) ERR;
   if (mode != expected_mode) ERR;
   if (extended_format != expected_extended_format) ERR;

   /* Nothing to do with inq_format, but let's check the base_pe
    * functions. */
   if (expected_format == NC_FORMAT_CLASSIC || expected_format == NC_FORMAT_64BIT_OFFSET ||
       expected_format == NC_FORMAT_CDF5) {
      if (nc_set_base_pe(ncid, 0)) ERR;
      if (nc_inq_base_pe(ncid, NULL)) ERR;
   } else {
      if (nc_set_base_pe(ncid, 0) != NC_ENOTNC3) ERR;
      if (nc_inq_base_pe(ncid, NULL) != NC_ENOTNC3) ERR;
   }

   return 0;
}

int
main(int argc, char **argv)
{
   printf("\n*** Testing HDF4/NetCDF-4 interoperability...\n");
   printf("*** testing that all hdf4 files can be opened...");
   {
#define NUM_SAMPLE_FILES 5
      int ncid;
      int nvars_in, ndims_in, natts_in, unlimdim_in;
      char file_name[NUM_SAMPLE_FILES][NC_MAX_NAME + 1] = {"AMSR_E_L2_Rain_V10_200905312326_A.hdf",
							   "AMSR_E_L3_DailyLand_V06_20020619.hdf",
							   "MOD29.A2000055.0005.005.2006267200024.hdf",
							   "MYD29.A2002185.0000.005.2007160150627.hdf",
							   "MYD29.A2009152.0000.005.2009153124331.hdf"};
      int expected_mode = NC_NETCDF4;
      int expected_extended_format = NC_FORMATX_NC_HDF4;
      int f;

      for (f = 0; f < NUM_SAMPLE_FILES; f++)
      {
	 if (nc_open(file_name[f], NC_NOWRITE, &ncid)) ERR;
	 if (nc_inq(ncid, &ndims_in, &nvars_in, &natts_in, &unlimdim_in)) ERR;
         if (check_inq_format(ncid, NC_FORMAT_NETCDF4, expected_extended_format, expected_mode)) ERR;
	 if (nc_close(ncid)) ERR;
      }
   }
   SUMMARIZE_ERR;
   FINAL_RESULTS;
}
