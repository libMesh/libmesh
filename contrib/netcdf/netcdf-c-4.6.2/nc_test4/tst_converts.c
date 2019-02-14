/* This is part of the netCDF package.
   Copyright 2005 University Corporation for Atmospheric Research/Unidata
   See COPYRIGHT file for conditions of use.

   Test data conversions and fill value handling.

   $Id: tst_converts.c,v 1.16 2008/10/20 01:48:09 ed Exp $
*/

#include <nc_tests.h>
#include "err_macros.h"
#include "netcdf.h"

#define FILE_NAME "tst_converts.nc"
#define ATT1_NAME "att1"
#define ATT2_NAME "att2"
#define DIM1_NAME "dim1"
#define DIM1_LEN 3
#define DIM2_NAME "dim2"
#define DIM2_LEN 15
#define VAR1_NAME "var1"
#define VAR2_NAME "var2"

/* This is handy for print statements. */
static char *format_name[MAX_NUM_FORMATS] = {"classic", "64-bit offset", "netCDF-4",
                                             "netCDF-4 classic model", "CDF5"};

int check_file(int format, unsigned char *uchar_out);
int create_file(int format, unsigned char *uchar_out);

/* Determine how many formats are available, and what they are. */
void
determine_test_formats(int *num_formats, int *format)
{
   int ind = 0;
   int num;

   /* Check inputs. */
   assert(num_formats && format);

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
   unsigned char uchar_out[DIM1_LEN] = {0, 128, 255};
   int format[MAX_NUM_FORMATS];
   int num_formats;
   int f = 0;

   printf("\n*** Testing netcdf data conversion.\n");
   determine_test_formats(&num_formats, format);

   for (f = 0; f < num_formats; f++)
   {
      printf("*** Testing conversion in netCDF %s files... ", format_name[f]);
      create_file(format[f], uchar_out);
      check_file(format[f], uchar_out);
      SUMMARIZE_ERR;
   }

   FINAL_RESULTS;
}

/* Create a test file with one var of type NC_BYTE. */
int
create_file(int format, unsigned char *uchar_out)
{
   int ncid, varid, cflags=0, dimids[1];
   int retval;

   if (format == NC_FORMAT_64BIT_OFFSET)
      cflags |= NC_64BIT_OFFSET;
   else if (format == NC_FORMAT_CDF5)
      cflags |= NC_CDF5;
   else if (format == NC_FORMAT_NETCDF4_CLASSIC)
   {
      cflags |= (NC_NETCDF4|NC_CLASSIC_MODEL);
   }
   else if (format == NC_FORMAT_NETCDF4)
      cflags |= NC_NETCDF4;

   if (nc_create(FILE_NAME, cflags, &ncid)) ERR;
   if (nc_def_dim(ncid, DIM1_NAME, DIM1_LEN, &dimids[0])) ERR;
   if (nc_def_var(ncid, VAR1_NAME, NC_BYTE, 1, dimids, &varid)) ERR;
   if (nc_enddef(ncid)) ERR;
   retval = nc_put_var_uchar(ncid, varid, uchar_out);
   if (format == NC_FORMAT_NETCDF4 || format == NC_FORMAT_64BIT_DATA)
   {
      if (retval != NC_ERANGE) ERR;
   }
   else if (retval != NC_NOERR) ERR;

   if (nc_close(ncid)) ERR;
   return NC_NOERR;
}

int
check_file(int format, unsigned char *uchar_out)
{
   int ncid;
   int ndims, natts;
   int dimids_var[1], var_type;
   char var_name[NC_MAX_NAME+1];
   unsigned char uchar_in[DIM1_LEN];
   signed char char_in[DIM1_LEN];
   unsigned short ushort_in[DIM1_LEN];
   short short_in[DIM1_LEN];
   unsigned int uint_in[DIM1_LEN];
   int int_in[DIM1_LEN];
   long long int64_in[DIM1_LEN];
   unsigned long long uint64_in[DIM1_LEN];
   int i, res;

   /* Read it back in, and check conversions. */
   if (nc_open(FILE_NAME, NC_NOWRITE, &ncid)) ERR;
   if (nc_inq_var(ncid, 0, var_name, &var_type, &ndims, dimids_var, &natts)) ERR;
   if (strcmp(var_name, VAR1_NAME) || natts !=0 || ndims != 1 ||
       dimids_var[0] != 0 || var_type != NC_BYTE) ERR;

   /* This is actually an NC_BYTE, with some negatives, so this should
    * generate a range error for netcdf-4, but not for netcdf-3,
    * because range errors are not generated for byte type
    * conversions. */
   res = nc_get_var_uchar(ncid, 0, uchar_in);
   if (format == NC_FORMAT_NETCDF4 || format == NC_FORMAT_64BIT_DATA)
   {
      if (res != NC_ERANGE) ERR;
   }
   else if (res) ERR;

   for (i=0; i<DIM1_LEN; i++)
#ifdef ERANGE_FILL
      if (uchar_in[i] != uchar_out[i] && uchar_in[i] != NC_FILL_UBYTE) ERR;
#else
      if (uchar_in[i] != uchar_out[i]) ERR;
#endif

   if (nc_get_var_schar(ncid, 0, char_in)) ERR;
   for (i=0; i<DIM1_LEN; i++)
#ifdef ERANGE_FILL
      if (char_in[i] != (signed char)uchar_out[i] && char_in[i] != NC_FILL_BYTE) ERR;
#else
      if (char_in[i] != (signed char)uchar_out[i]) ERR;
#endif

   if (nc_get_var_short(ncid, 0, short_in)) ERR;
   for (i=0; i<DIM1_LEN; i++)
#ifdef ERANGE_FILL
      if (short_in[i] != (signed char)uchar_out[i] && short_in[i] != NC_FILL_BYTE) ERR;
#else
      if (short_in[i] != (signed char)uchar_out[i]) ERR;
#endif

   if (nc_get_var_int(ncid, 0, int_in)) ERR;
   for (i=0; i<DIM1_LEN; i++)
#ifdef ERANGE_FILL
      if (int_in[i] != (signed char)uchar_out[i] && int_in[i] != NC_FILL_BYTE) ERR;
#else
      if (int_in[i] != (signed char)uchar_out[i]) ERR;
#endif

   if (format == NC_FORMAT_NETCDF4 || format == NC_FORMAT_NETCDF4_CLASSIC)
   {
      /* Since we wrote them as NC_BYTE, some of these are negative
       * values, and will return a range error when reading into
       * unsigned type. To compare values, first cast uchar_out to
       * signed int, then cast again to the type we are reading it
       * as. */
      if (nc_get_var_ushort(ncid, 0, ushort_in) != NC_ERANGE) ERR;
      for (i=0; i<DIM1_LEN; i++)
	 if (ushort_in[i] != (unsigned short)(signed char)uchar_out[i]) ERR;

      if (nc_get_var_uint(ncid, 0, uint_in) != NC_ERANGE) ERR;
      for (i=0; i<DIM1_LEN; i++)
	 if (uint_in[i] != (unsigned int)(signed char)uchar_out[i]) ERR;

      if (nc_get_var_ulonglong(ncid, 0, uint64_in) != NC_ERANGE) ERR;
      for (i=0; i<DIM1_LEN; i++)
	 if (uint64_in[i] != (unsigned long long)(signed char)uchar_out[i]) ERR;

      if (nc_get_var_longlong(ncid, 0, int64_in)) ERR;
      for (i=0; i<DIM1_LEN; i++)
	 if (int64_in[i] != (signed char)uchar_out[i]) ERR;

   }

   if (nc_close(ncid)) ERR;
   return 0;
}
