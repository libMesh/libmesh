/* This is part of the netCDF package.
   Copyright 2005 University Corporation for Atmospheric Research/Unidata
   See COPYRIGHT file for conditions of use.

   Test netcdf-4 types.

   Ed Hartnett
*/

#include "config.h"
#include "nc_tests.h"
#include "err_macros.h"

#define MAX_VARNAME 20
#define NUM_TYPES 6
#define NUM_DIMS 1
#define SIZE 5
#define STRIDE_SIZE 2
#define FILENAME "tst_types.nc"
#define FILENAME2 "tst_types2.nc"
#define FILENAME3 "tst_types3.nc"
#define FILENAME4 "tst_types4.nc"
#define A_NAME "some_name"

#define CLEAN_INPUT_BUFFERS                     \
   for (i = 0; i < SIZE; i++) {                 \
      ubyte_data_out[i] = 0;                    \
      ushort_data_out[i] = 0;                   \
      uint_data_out[i] = 0;                     \
      int64_data_out[i] = 0;                    \
      uint64_data_out[i] = 0;                   \
   }

/* Create the test file for these tests. */
int
create_test_file(char *filename, int *varid, int *ncid)
{
   int dimid;
   int type;
   char varname[MAX_VARNAME];

   /* Open a netcdf-4 file, and one dimension. */
   if (nc_create(filename, NC_NETCDF4, ncid)) ERR;
   if (nc_def_dim(*ncid, "dim1", SIZE, &dimid)) ERR;

   /* Create vars of the new types. Take advantage of the fact that
    * new types are numbered from NC_UBYTE (7) through NC_STRING
    * (12). */
   for (type = 0; type < NUM_TYPES; type++)
   {
      /* Create a var... */
      sprintf(varname, "var_%d", type);
      if (nc_def_var(*ncid, varname, type + NC_UBYTE, 1, &dimid, &varid[type])) ERR;
   }
   return 0;
}

int main(int argc, char *argv[])
{
   /* IDs, names, and parameters for the var[asm1] functions. */
   int ncid, varid[NUM_TYPES], dimid;
   char varname[MAX_VARNAME];
   size_t index1[NUM_DIMS], start[NUM_DIMS];
   size_t count[NUM_DIMS];
   ptrdiff_t imap[NUM_DIMS];
   ptrdiff_t stride[NUM_DIMS];

   /* Phoney data we will write. */
   unsigned char ubyte_data_out[] = {0,1,2,3,4};
   unsigned short ushort_data_out[] = {0,11,22,33,44};
   unsigned int uint_data_out[] = {0,111,222,333,3000000000u};
   long long int int64_data_out[] = {0,-111111111,2222222222,-3333333333,444444444};
   unsigned long long int uint64_data_out[] = {0,111111111,2222222222,33333333,44444444};

   /* We will read back in the phoney data with these. */
   unsigned char ubyte_data_in[SIZE];
   unsigned short ushort_data_in[SIZE];
   unsigned int uint_data_in[SIZE];
   long long int int64_data_in[SIZE];
   unsigned long long int uint64_data_in[SIZE];

   int i;
   int type;

   printf("\n*** Testing netCDF-4 types...\n");
   printf("*** testing netCDF-4 new atomic types and equality...");
   {
      int ncid1, ncid2;
      int equal;

      /* Open a netcdf-4 file, and one dimension. */
      if (nc_create(FILENAME, NC_NETCDF4, &ncid1)) ERR;
      if (nc_def_dim(ncid1, "dim1", SIZE, &dimid)) ERR;

      /* Open another netcdf-4 file. */
      if (nc_create(FILENAME2, NC_NETCDF4, &ncid2)) ERR;

      /* Create vars of the new types. Take advantage of the fact that
       * new types are numbered from NC_UBYTE (7) through NC_STRING
       * (12). */
      for (type = 0; type < NUM_TYPES; type++)
      {
         /* Create a var... */
         sprintf(varname, "var_%d", type);
         if (nc_def_var(ncid1, varname, type + NC_UBYTE, 1, &dimid, &varid[type])) ERR;
         if (nc_inq_type_equal(ncid1, type + NC_UBYTE, ncid2, type + NC_UBYTE, &equal)) ERR;
         if (!equal) ERR;
      }
      nc_close(ncid2);
      nc_close(ncid1);
   }
   SUMMARIZE_ERR;

#define ENUM_TYPE_NAME "enum_type"
#define ENUM_FIELD1_NAME "enum_field_1"
#define ENUM_UNEQUAL_TYPE_NAME_3 "enum_type_3"
#define ENUM_UNEQUAL_TYPE_NAME_4 "enum_type_4"
#define NUM_CLASSIC_TYPES 6
#define NUM_ENHANCED_TYPES 12

   printf("*** testing user-defined netCDF-4 enum types and equality...");
   {
      int ncid1, ncid2, ncid3, ncid4;
      int typeid1, typeid2, typeid3;
      int enum_value = TEST_VAL_42;
      int equal;
      int classic_type[NUM_CLASSIC_TYPES] = {NC_BYTE, NC_CHAR, NC_SHORT, NC_INT, NC_FLOAT, NC_DOUBLE};
      int enhanced_type[NUM_ENHANCED_TYPES] = {NC_BYTE, NC_CHAR, NC_SHORT, NC_INT, NC_FLOAT, NC_DOUBLE,
                                               NC_UBYTE, NC_USHORT, NC_UINT, NC_INT64, NC_UINT64, NC_STRING};
      int t;

      /* Create a netcdf-4 file. */
      if (nc_create(FILENAME, NC_NETCDF4, &ncid1)) ERR;

      /* Create an enum type. */
      if (nc_def_enum(ncid1, NC_INT, ENUM_TYPE_NAME, &typeid1)) ERR;
      if (nc_insert_enum(ncid1, typeid1, ENUM_FIELD1_NAME, &enum_value)) ERR;

      /* Create another netcdf-4 file. */
      if (nc_create(FILENAME2, NC_NETCDF4, &ncid2)) ERR;

      /* Create a netcdf-3 classic file. */
      if (nc_create(FILENAME3, 0, &ncid3)) ERR;

      /* Create a netcdf-4 classic model file. */
      if (nc_create(FILENAME4, NC_NETCDF4|NC_CLASSIC_MODEL, &ncid4)) ERR;

      /* Create an enum type that will be equal to typeid1. */
      if (nc_def_enum(ncid2, NC_INT, ENUM_TYPE_NAME, &typeid2)) ERR;
      if (nc_insert_enum(ncid2, typeid2, ENUM_FIELD1_NAME, &enum_value)) ERR;

      /* These will fail. */
      if (nc_def_enum(ncid2, NC_INT, ENUM_TYPE_NAME, &typeid3) != NC_ENAMEINUSE) ERR;
      if (nc_def_enum(ncid2 + TEST_VAL_42, NC_INT, ENUM_UNEQUAL_TYPE_NAME_3, &typeid3) != NC_EBADID) ERR;
      if (nc_def_enum(ncid2, NC_INT, NULL, &typeid3) != NC_EINVAL) ERR;
      if (nc_def_enum(ncid3, NC_SHORT, ENUM_UNEQUAL_TYPE_NAME_3, &typeid3) != NC_ENOTNC4) ERR;
      if (nc_def_enum(ncid4, NC_SHORT, ENUM_UNEQUAL_TYPE_NAME_3, &typeid3) != NC_ESTRICTNC3) ERR;
      if (nc_def_opaque(ncid4, TEST_VAL_42, A_NAME, &typeid3) != NC_ESTRICTNC3) ERR;
      if (nc_def_compound(ncid4, TEST_VAL_42, A_NAME, &typeid3) != NC_ESTRICTNC3) ERR;
      if (nc_def_vlen(ncid4, A_NAME, NC_INT, &typeid3) != NC_ESTRICTNC3) ERR;

      /* Create some enum types that will not be equal to typeid1. */
      if (nc_def_enum(ncid2, NC_SHORT, ENUM_UNEQUAL_TYPE_NAME_3, &typeid3)) ERR;
      if (nc_insert_enum(ncid2, typeid3, ENUM_FIELD1_NAME, &enum_value)) ERR;

      /* This will fail because types are not yet committed. */
      if (nc_inq_type_equal(ncid1, typeid1, ncid2, typeid2, &equal) != NC_EHDFERR) ERR;

      /* Commit the types to the file. */
      if (nc_enddef(ncid1)) ERR;
      if (nc_enddef(ncid2)) ERR;

      /* This succeeds but does nothing. */
      if (nc_inq_type_equal(ncid1, typeid1, ncid2, typeid2, NULL)) ERR;

      /* These will fail. */
      if (nc_inq_type_equal(ncid1, 0, ncid2, typeid2, &equal) != NC_EINVAL) ERR;
      if (nc_inq_type_equal(ncid1, typeid1, ncid2, 0, &equal) != NC_EINVAL) ERR;
      if (nc_inq_type_equal(ncid1, -TEST_VAL_42, ncid2, typeid2, &equal) != NC_EINVAL) ERR;
      if (nc_inq_type_equal(ncid1, typeid1, ncid2, -TEST_VAL_42, &equal) != NC_EINVAL) ERR;
      if (nc_inq_type_equal(ncid1 + TEST_VAL_42, typeid1, ncid2, typeid2, &equal) != NC_EBADID) ERR;
      if (nc_inq_type_equal(ncid1, typeid1 + TEST_VAL_42, ncid2, typeid2, &equal) != NC_EBADTYPE) ERR;
      if (nc_inq_type_equal(ncid1, typeid1, ncid2 + TEST_VAL_42, typeid2, &equal) != NC_EBADID) ERR;
      if (nc_inq_type_equal(ncid1, typeid1, ncid2, typeid2 + TEST_VAL_42, &equal) != NC_EBADTYPE) ERR;

      /* Ensure the two equal types are equal. */
      if (nc_inq_type_equal(ncid1, typeid1, ncid2, typeid2, &equal)) ERR;
      if (!equal) ERR;

      /* Ensure the two unequal types are not equal. */
      if (nc_inq_type_equal(ncid1, typeid1, ncid2, typeid3, &equal)) ERR;
      if (equal) ERR;

      /* Atomic types are not equal to user-defined types, but are
       * equal to themselves. */
      for (t = 0; t < NUM_CLASSIC_TYPES; t++)
      {
         if (nc_inq_type_equal(ncid1, typeid1, ncid3, classic_type[t], &equal)) ERR;
         if (equal) ERR;
         if (nc_inq_type_equal(ncid1, classic_type[t], ncid3, classic_type[t], &equal)) ERR;
         if (!equal) ERR;
      }
      
      for (t = 0; t < NUM_ENHANCED_TYPES; t++)
      {
         if (nc_inq_type_equal(ncid1, typeid1, ncid2, enhanced_type[t], &equal)) ERR;
         if (equal) ERR;
         if (nc_inq_type_equal(ncid1, enhanced_type[t], ncid2, enhanced_type[t], &equal)) ERR;
         if (!equal) ERR;
      }

      /* Close the files. */
      nc_close(ncid1);
      nc_close(ncid2);
      nc_close(ncid3);
      nc_close(ncid4);
   }
   SUMMARIZE_ERR;

   /* Test the vara functions. */
   printf("*** testing vara functions with new netCDF-4 atomic types...");
   {
      CLEAN_INPUT_BUFFERS;
      start[0] = 0;
      count[0] = SIZE;

      /* Open a netcdf-4 file, and one dimension. */
      if (create_test_file(FILENAME, varid, &ncid)) ERR;

      /* This will not work. */
      if (nc_put_vara_uchar(ncid, varid[0], NULL, count, ubyte_data_out) !=
          NC_EINVALCOORDS) ERR;

      /* This will work. */
      if (nc_put_vara_uchar(ncid, varid[0], start, count, ubyte_data_out)) ERR;

      /* This will not work. */
      if (nc_get_vara_uchar(ncid, varid[0], NULL, count, ubyte_data_in) !=
          NC_EINVALCOORDS) ERR;

      /* This will work. */
      if (nc_get_vara_uchar(ncid, varid[0], start, NULL, ubyte_data_in)) ERR;
      for (i = 0; i < SIZE; i++)
         if (ubyte_data_in[i] != ubyte_data_out[i]) ERR;

      CLEAN_INPUT_BUFFERS;
      if (nc_get_vara_uchar(ncid, varid[0], start, count, ubyte_data_in)) ERR;
      for (i = 0; i < SIZE; i++)
         if (ubyte_data_in[i] != ubyte_data_out[i]) ERR;

      /* NULL count will be interpreted as count of full extent. */
      if (nc_put_vara_ushort(ncid, varid[1], start, NULL, ushort_data_out)) ERR;
      if (nc_get_vara_ushort(ncid, varid[1], start, count, ushort_data_in)) ERR;
      for (i = 0; i < SIZE; i++)
         if (ubyte_data_in[i] != ubyte_data_out[i]) ERR;

      /* vars with NULL stride is the same as vara. */
      if (nc_put_vars_uint(ncid, varid[2], start, NULL, NULL, uint_data_out)) ERR;
      if (nc_get_vara_uint(ncid, varid[2], start, count, uint_data_in)) ERR;
      for (i = 0; i < SIZE; i++)
         if (ubyte_data_in[i] != ubyte_data_out[i]) ERR;

      if (nc_put_vara_longlong(ncid, varid[3], start, count, int64_data_out)) ERR;
      if (nc_get_vars_longlong(ncid, varid[3], start, NULL, NULL, int64_data_in)) ERR;
      for (i = 0; i < SIZE; i++)
         if (ubyte_data_in[i] != ubyte_data_out[i]) ERR;

      if (nc_put_vara_ulonglong(ncid, varid[4], start, count, uint64_data_out)) ERR;
      if (nc_get_vara_ulonglong(ncid, varid[4], start, count, uint64_data_in)) ERR;
      for (i = 0; i < SIZE; i++)
         if (ubyte_data_in[i] != ubyte_data_out[i]) ERR;

      /* Close the test file. */
      nc_close(ncid);
   }
   SUMMARIZE_ERR;

   /* Test the vars functions. */
   printf("*** testing vars functions with new netCDF-4 atomic types...");
   {
      CLEAN_INPUT_BUFFERS;
      start[0] = 0;
      count[0] = 2;
      stride[0] = STRIDE_SIZE;

      /* Open a netcdf-4 file, and one dimension. */
      if (create_test_file(FILENAME, varid, &ncid)) ERR;

      if (nc_put_vars_uchar(ncid, varid[0], start, count, stride,
                            ubyte_data_out)) ERR;

      /* This will not work. */
      if (nc_get_vars_uchar(ncid, varid[0], NULL, count, stride,
                            ubyte_data_in) != NC_EINVALCOORDS) ERR;

      /* This will work. */
      if (nc_get_vars_uchar(ncid, varid[0], start, count, stride,
                            ubyte_data_in)) ERR;
      if (ubyte_data_in[0] != ubyte_data_out[0]) ERR;
      if (ubyte_data_in[1] != ubyte_data_out[STRIDE_SIZE]) ERR;

      /* This will not work. */
      if (nc_put_vars_ushort(ncid, varid[1], NULL, count, stride,
                             ushort_data_out) != NC_EINVALCOORDS) ERR;

      /* This will work. */
      if (nc_put_vars_ushort(ncid, varid[1], start, count, stride,
                             ushort_data_out)) ERR;
      if (nc_get_vars_ushort(ncid, varid[1], start, count, stride,
                             ushort_data_in)) ERR;
      for (i = 0; i < 2; i++)
         if (ubyte_data_in[i] != ubyte_data_out[i]) ERR;

      if (nc_put_vars_uint(ncid, varid[2], start, count, stride,
                           uint_data_out)) ERR;
      if (nc_get_vars_uint(ncid, varid[2], start, count, stride,
                           uint_data_in)) ERR;
      for (i = 0; i < 2; i++)
         if (ubyte_data_in[i] != ubyte_data_out[i]) ERR;

      if (nc_put_vars_longlong(ncid, varid[3], start, count, stride,
                               int64_data_out)) ERR;
      if (nc_get_vars_longlong(ncid, varid[3], start, count, stride,
                               int64_data_in)) ERR;
      for (i = 0; i < 2; i++)
         if (ubyte_data_in[i] != ubyte_data_out[i]) ERR;

      if (nc_put_vars_ulonglong(ncid, varid[4], start, count, stride,
                                uint64_data_out)) ERR;
      if (nc_get_vars_ulonglong(ncid, varid[4], start, count, stride,
                                uint64_data_in)) ERR;
      for (i = 0; i < 2; i++)
         if (ubyte_data_in[i] != ubyte_data_out[i]) ERR;

      /* Close the test file. */
      nc_close(ncid);
   }
   SUMMARIZE_ERR;

   /* Test the var1 functions. */
   printf("*** testing var1 functions with new netCDF-4 atomic types...");
   {
      CLEAN_INPUT_BUFFERS;
      index1[0] = 0;

      /* Open a netcdf-4 file, and one dimension. */
      if (create_test_file(FILENAME, varid, &ncid)) ERR;

      /* This will not work. */
      if (nc_put_var1_uchar(ncid, varid[0], NULL, ubyte_data_out) !=
          NC_EINVALCOORDS) ERR;

      /* This will work. */
      if (nc_put_var1_uchar(ncid, varid[0], index1, ubyte_data_out)) ERR;

      /* This will not work. */
      if (nc_get_var1_uchar(ncid, varid[0], NULL, ubyte_data_in) != NC_EINVALCOORDS) ERR;

      /* This will work. */
      if (nc_get_var1_uchar(ncid, varid[0], index1, ubyte_data_in)) ERR;
      if (ubyte_data_in[0] != ubyte_data_out[0]) ERR;

      if (nc_put_var1_ushort(ncid, varid[1], index1, ushort_data_out)) ERR;
      if (nc_get_var1_ushort(ncid, varid[1], index1, ushort_data_in)) ERR;
      if (ushort_data_in[0] != ushort_data_out[0]) ERR;

      if (nc_put_var1_uint(ncid, varid[2], index1, uint_data_out)) ERR;
      if (nc_get_var1_uint(ncid, varid[2], index1, uint_data_in)) ERR;
      if (uint_data_in[0] != uint_data_out[0]) ERR;

      if (nc_put_var1_longlong(ncid, varid[3], index1, int64_data_out)) ERR;
      if (nc_get_var1_longlong(ncid, varid[3], index1, int64_data_in)) ERR;
      if (int64_data_in[0] != int64_data_out[0]) ERR;

      if (nc_put_var1_ulonglong(ncid, varid[4], index1, uint64_data_out)) ERR;
      if (nc_get_var1_ulonglong(ncid, varid[4], index1, uint64_data_in)) ERR;
      if (uint64_data_in[0] != uint64_data_out[0]) ERR;

      /* Close the test file. */
      nc_close(ncid);
   }
   SUMMARIZE_ERR;

   /* Test the varm functions. */
   printf("*** testing varm functions with new netCDF-4 atomic types...");
   {
      int ncid, varid[NUM_TYPES];
      CLEAN_INPUT_BUFFERS;
      start[0] = 0;
      count[0] = 1;
      stride[0] = 1;
      imap[0] = 0;

      /* Open a netcdf-4 file, and one dimension. */
      if (create_test_file(FILENAME, varid, &ncid)) ERR;

      /* This will not work. */
      if (nc_put_varm_ubyte(ncid, varid[0], NULL, count, stride, imap,
                            ubyte_data_out) != NC_EINVALCOORDS) ERR;

      /* This will work. */
      if (nc_put_varm_ubyte(ncid, varid[0], start, count, stride, imap,
                            ubyte_data_out)) ERR;

      /* This will not work. */
      if (nc_get_varm_ubyte(ncid, varid[0], NULL, count, stride, imap,
                            ubyte_data_in) != NC_EINVALCOORDS) ERR;

      /* This will work. */
      if (nc_get_varm_ubyte(ncid, varid[0], start, count, stride, imap,
                            ubyte_data_in)) ERR;
      for (i = 0; i < STRIDE_SIZE; i++)
         if (ubyte_data_in[i] != ubyte_data_out[i]) ERR;

      if (nc_put_varm_ushort(ncid, varid[1], start, count, stride, imap,
                             ushort_data_out)) ERR;
      if (nc_get_varm_ushort(ncid, varid[1], start, count, stride, imap,
                             ushort_data_in)) ERR;
      for (i = 0; i < STRIDE_SIZE; i++)
         if (ubyte_data_in[i] != ubyte_data_out[i]) ERR;

      if (nc_put_varm_uint(ncid, varid[2], start, count, stride, imap,
                           uint_data_out)) ERR;
      if (nc_get_varm_uint(ncid, varid[2], start, count, stride, imap,
                           uint_data_in)) ERR;
      for (i = 0; i < STRIDE_SIZE; i++)
         if (ubyte_data_in[i] != ubyte_data_out[i]) ERR;

      if (nc_put_varm_longlong(ncid, varid[3], start, count,
                               stride, imap, int64_data_out)) ERR;
      if (nc_get_varm_longlong(ncid, varid[3], start, count,
                               stride, imap, int64_data_in)) ERR;
      for (i = 0; i < STRIDE_SIZE; i++)
         if (ubyte_data_in[i] != ubyte_data_out[i]) ERR;

      if (nc_put_varm_ulonglong(ncid, varid[4], start, count,
                                stride, imap, uint64_data_out)) ERR;
      if (nc_get_varm_ulonglong(ncid, varid[4], start, count,
                                stride, imap, uint64_data_in)) ERR;
      for (i = 0; i < STRIDE_SIZE; i++)
         if (ubyte_data_in[i] != ubyte_data_out[i]) ERR;

      /* Close the test file. */
      nc_close(ncid);
   }
   SUMMARIZE_ERR;

   FINAL_RESULTS;
}
