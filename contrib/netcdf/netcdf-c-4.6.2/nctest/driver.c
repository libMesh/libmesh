/*********************************************************************
 * Copyright 1993-2006, UCAR/Unidata
 * See COPYRIGHT file for copying and redistribution conditions.
 *
 * Test driver for netCDF implementation.  This program performs tests
 * against the netCDF specification for all user-level functions in an
 * implementation of the netCDF library.  Must be invoked from a
 * directory in which the invoker has write permission.
 *
 * Glenn Davis, Russ Rew, Ed Hartnett
 *********************************************************************/

#include <config.h>
#include <stdlib.h>
#include <netcdf.h>
#include <stdio.h>
#include <string.h>
#include "testcdf.h"		/* defines in-memory test netcdf structure */
#include "tests.h"

#define MAX_NUM_FORMATS 5

int current_format = 0;

/* Determine how many formats are available, and what they are. */
void
determine_test_formats(int *num_formats, int *format)
{
   int ind = 0;
   int num;

   /* We always have classic and 64-bit offset */
   num = 2;
   format[ind++] = NC_FORMAT_CLASSIC;
   format[ind++] = NC_FORMAT_64BIT_OFFSET;

   /* Do we have netCDF-4 and netCDF-4 classic? */
#ifdef USE_NETCDF4
   num += 2;
   format[ind++] = NC_FORMAT_NETCDF4_CLASSIC;
   format[ind++] = NC_FORMAT_NETCDF4;
#endif /* USE_NETCDF4 */

   /* Do we have CDF5? */
#ifdef ENABLE_CDF5
   num++;
   format[ind++] = NC_FORMAT_CDF5;
#endif /* ENABLE_CDF5 */

   *num_formats = num;
}

int
main(int argc, char **argv)
{
   /*EXTERNL int ncopts;	*/	/* netCDF error options */
   char *format_name[MAX_NUM_FORMATS] = {"classic", "64bit_offset", "netcdf4_classic",
                                         "netcdf4", "CDF5"};
   char testfile[NC_MAX_NAME];
   int format[MAX_NUM_FORMATS];
   int num_formats;
   int i, nerrs = 0;

   ncopts &= ~NC_FATAL;	/* make errors nonfatal */
   ncopts &= ~NC_VERBOSE;	/* turn off error messages */

   /* How many formats are we testing? */
   determine_test_formats(&num_formats, format);
   printf("Testing V2 API with %d different netCDF formats.\n", num_formats);

   for (i = 0; i < num_formats; i++)
   {
      current_format = format[i];

      /* Skip netCDF-4 - only netCDF-4 classic will work. */
      if (format[i] == NC_FORMAT_NETCDF4)
         continue;
      
      /* Come up with a test file name. */
      sprintf(testfile, "nctest_%s.nc", format_name[i]);
      printf("Testing %s with file %s.\n", format_name[i], testfile);      

      /* Set the default format. */
      nc_set_default_format(format[i], NULL);

      /* Run all the tests for this format. */
      nerrs += test_nccreate(testfile);
      nerrs += test_ncopen(testfile);
      nerrs += test_ncredef(testfile);
      nerrs += test_ncendef(testfile);
      nerrs += test_ncclose(testfile);
      nerrs += test_ncinquire(testfile);
      nerrs += test_ncsync(testfile);
      nerrs += test_ncabort(testfile);
      nerrs += test_ncdimdef(testfile);
      nerrs += test_ncdimid(testfile);
      nerrs += test_ncdiminq(testfile);
      nerrs += test_ncdimrename(testfile);
      nerrs += test_ncvardef(testfile);
      nerrs += test_ncvarid(testfile);
      nerrs += test_ncvarinq(testfile);
      nerrs += test_ncvarputg(testfile);
      nerrs += test_ncvarput1(testfile);
      nerrs += test_ncvarget1(testfile);
      nerrs += test_ncvarput(testfile);
      nerrs += test_ncvarget(testfile);
      nerrs += test_ncvarputg(testfile);
      nerrs += test_ncvargetg(testfile);
      nerrs += test_ncrecinq(testfile);
      nerrs += test_ncrecput(testfile);
      nerrs += test_ncrecget(testfile);
      nerrs += test_ncvarrename(testfile);
      nerrs += test_ncattput(testfile);
      nerrs += test_ncattinq(testfile);
      nerrs += test_ncattget(testfile);
      nerrs += test_ncattcopy(testfile, "test2.nc");
      nerrs += test_ncattname(testfile);
      nerrs += test_ncattrename(testfile);
      nerrs += test_ncattdel(testfile);
      nerrs += test_nctypelen();

      /* Clean up in-memory struct. */
      {
         int i;

         for (i = 0; i < test.ndims; i++)
            free(test.dims[i].name);

         for (i = 0; i < test.nvars; i++)
         {
            free(test.vars[i].name);
            free(test.vars[i].dims);
         }

         for (i = 0; i < test.natts; i++)
            free(test.atts[i].name);

      }
   }
    
   fprintf(stderr, "\nTotal number of failures: %d\n", nerrs);

   if (nerrs)
   {
      fprintf(stderr, "nctest FAILURE!!!\n");
      return 2;
   }
   else
      fprintf(stderr, "nctest SUCCESS!!!\n");

   return 0;
}
