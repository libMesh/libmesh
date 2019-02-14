/* Copyright 2007-2011, UCAR/Unidata. See COPYRIGHT file for copying
 * and redistribution conditions.
 *
 * This is part of the netCDF package.
 *
 * This test is for parallel IO and the collective access of metadata
 * with HDF5.
 *
 * Ward Fisher, Ed Hartnett
 */

#include "config.h"
#include "nc_tests.h"
#include "err_macros.h"

#define TEST_NAME "tst_parallel4_simplerw_coll"
#define NDIMS 3
#define DIMSIZE 16
#define NUM_SLABS 16
#define DIM1_NAME "slab"
#define DIM2_NAME "x"
#define DIM3_NAME "y"
#define VAR_NAME "Bond_James_Bond"
#define NUM_FILL_TEST_RUNS 3

int
main(int argc, char **argv)
{
   int mpi_namelen;
   char mpi_name[MPI_MAX_PROCESSOR_NAME];
   int mpi_size, mpi_rank;
   MPI_Comm comm = MPI_COMM_WORLD;
   MPI_Info info = MPI_INFO_NULL;
   double start_time = 0, total_time;
   int mpi_size_in;
#define NUM_TEST_TYPES 11
   nc_type test_type[NUM_TEST_TYPES] = {NC_BYTE, NC_CHAR, NC_SHORT, NC_INT, NC_FLOAT, NC_DOUBLE,
                                        NC_UBYTE, NC_USHORT, NC_UINT, NC_INT64, NC_UINT64};
   int tt, fv;
   int j, i, k, ret;

   /* Initialize MPI. */
   MPI_Init(&argc,&argv);
   MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
   MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
   MPI_Get_processor_name(mpi_name, &mpi_namelen);

   /* Must be able to evenly divide my slabs between processors. */
   if (NUM_SLABS % mpi_size)
   {
      if (!mpi_rank)
         printf("NUM_SLABS (%d) is not evenly divisible by mpi_size(%d)\n",
                NUM_SLABS, mpi_size);
      ERR;
   }

   if (!mpi_rank)
      printf("\n*** Testing parallel I/O some more.\n");

   /* Test for different fill value settings. */
   for (fv = 0; fv < NUM_FILL_TEST_RUNS; fv++)
   {
      /* Test for different netCDF types. */
      for (tt = 0; tt < NUM_TEST_TYPES; tt++)
      {
         char file_name[NC_MAX_NAME + 1];
         int fill_mode_in;
         void *data, *data_in;
         void *fill_value, *fill_value_in;
         size_t type_size;
         size_t write_start[NDIMS] = {0, 0, 1};
         size_t write_count[NDIMS] = {1, DIMSIZE, DIMSIZE - 1};
         size_t read_start[NDIMS] = {0, 0, 0};
         size_t read_count[NDIMS] = {1, DIMSIZE, DIMSIZE};
         int ncid, varid, dimids[NDIMS];
         int ndims_in, nvars_in, natts_in, unlimdimid_in;

         /* Fill values to be expected. */
         signed char byte_expected_fill_value;
         unsigned char char_expected_fill_value;
         short short_expected_fill_value;
         int int_expected_fill_value;
         float float_expected_fill_value;
         double double_expected_fill_value;
         unsigned char ubyte_expected_fill_value;
         unsigned short ushort_expected_fill_value;
         unsigned int uint_expected_fill_value;
         long long int int64_expected_fill_value;
         unsigned long long int uint64_expected_fill_value;

         /* Fill values used when writing. */
         signed char byte_fill_value = -TEST_VAL_42;
         unsigned char char_fill_value = 'x';
         short short_fill_value = TEST_VAL_42 * 100;
         int int_fill_value = TEST_VAL_42 * 1000;
         float float_fill_value = TEST_VAL_42 * 1000;
         double double_fill_value = TEST_VAL_42 * 1000;
         unsigned char ubyte_fill_value = TEST_VAL_42;
         unsigned short ushort_fill_value = TEST_VAL_42 * 100;
         unsigned int uint_fill_value = TEST_VAL_42 * 1000;
         long long int int64_fill_value = TEST_VAL_42 * 1000;
         unsigned long long int uint64_fill_value = TEST_VAL_42 * 1000;

         /* Fill values read in. */
         signed char byte_fill_value_in;
         unsigned char char_fill_value_in;
         short short_fill_value_in;
         int int_fill_value_in;
         float float_fill_value_in;
         double double_fill_value_in;
         unsigned char ubyte_fill_value_in;
         unsigned short ushort_fill_value_in;
         unsigned int uint_fill_value_in;
         long long int int64_fill_value_in;
         unsigned long long int uint64_fill_value_in;
         
         /* Data to write and read. */
         signed char byte_data[DIMSIZE * DIMSIZE], byte_data_in[DIMSIZE * DIMSIZE];
         unsigned char char_data[DIMSIZE * DIMSIZE], char_data_in[DIMSIZE * DIMSIZE];
         short short_data[DIMSIZE * DIMSIZE], short_data_in[DIMSIZE * DIMSIZE];
         int int_data[DIMSIZE * DIMSIZE], int_data_in[DIMSIZE * DIMSIZE];
         float float_data[DIMSIZE * DIMSIZE], float_data_in[DIMSIZE * DIMSIZE];
         double double_data[DIMSIZE * DIMSIZE], double_data_in[DIMSIZE * DIMSIZE];
         unsigned char ubyte_data[DIMSIZE * DIMSIZE], ubyte_data_in[DIMSIZE * DIMSIZE];
         unsigned short ushort_data[DIMSIZE * DIMSIZE], ushort_data_in[DIMSIZE * DIMSIZE];
         unsigned int uint_data[DIMSIZE * DIMSIZE], uint_data_in[DIMSIZE * DIMSIZE];
         long long int int64_data[DIMSIZE * DIMSIZE], int64_data_in[DIMSIZE * DIMSIZE];
         unsigned long long int uint64_data[DIMSIZE * DIMSIZE], uint64_data_in[DIMSIZE * DIMSIZE];
         
         if (!mpi_rank)
            printf("*** writing a %d x %d x %d file from %d processors for fill value test %d type %d...\n",
                   NUM_SLABS, DIMSIZE, DIMSIZE, mpi_size, fv, test_type[tt]);

         /* Initialize test data. */
         switch(test_type[tt])
         {
         case NC_BYTE:
            for (i = 0; i < DIMSIZE * DIMSIZE; i++)
               byte_data[i] = mpi_rank;
            data = byte_data;
            data_in = byte_data_in;
            byte_expected_fill_value = fv ? byte_fill_value : NC_FILL_BYTE;
            fill_value = &byte_expected_fill_value;
            fill_value_in = &byte_fill_value_in;
            break;
         case NC_CHAR:
            for (i = 0; i < DIMSIZE * DIMSIZE; i++)
               char_data[i] = mpi_rank;
            data = char_data;
            data_in = char_data_in;
            char_expected_fill_value = fv ? char_fill_value : NC_FILL_CHAR;
            fill_value = &char_expected_fill_value;
            fill_value_in = &char_fill_value_in;
            break;
         case NC_SHORT:
            for (i = 0; i < DIMSIZE * DIMSIZE; i++)
               short_data[i] = mpi_rank;
            data = short_data;
            data_in = short_data_in;
            short_expected_fill_value = fv ? short_fill_value : NC_FILL_SHORT;
            fill_value = &short_expected_fill_value;
            fill_value_in = &short_fill_value_in;
            break;
         case NC_INT:
            for (i = 0; i < DIMSIZE * DIMSIZE; i++)
               int_data[i] = mpi_rank;
            data = int_data;
            data_in = int_data_in;
            int_expected_fill_value = fv ? int_fill_value : NC_FILL_INT;
            fill_value = &int_expected_fill_value;
            fill_value_in = &int_fill_value_in;
            break;
         case NC_FLOAT:
            for (i = 0; i < DIMSIZE * DIMSIZE; i++)
               float_data[i] = mpi_rank;
            data = float_data;
            data_in = float_data_in;
            float_expected_fill_value = fv ? float_fill_value : NC_FILL_FLOAT;
            fill_value = &float_expected_fill_value;
            fill_value_in = &float_fill_value_in;
            break;
         case NC_DOUBLE:
            for (i = 0; i < DIMSIZE * DIMSIZE; i++)
               double_data[i] = mpi_rank;
            data = double_data;
            data_in = double_data_in;
            double_expected_fill_value = fv ? double_fill_value : NC_FILL_DOUBLE;
            fill_value = &double_expected_fill_value;
            fill_value_in = &double_fill_value_in;
            break;
         case NC_UBYTE:
            for (i = 0; i < DIMSIZE * DIMSIZE; i++)
               ubyte_data[i] = mpi_rank;
            data = ubyte_data;
            data_in = ubyte_data_in;
            ubyte_expected_fill_value = fv ? ubyte_fill_value : NC_FILL_UBYTE;
            fill_value = &ubyte_expected_fill_value;
            fill_value_in = &ubyte_fill_value_in;
            break;
         case NC_USHORT:
            for (i = 0; i < DIMSIZE * DIMSIZE; i++)
               ushort_data[i] = mpi_rank;
            data = ushort_data;
            data_in = ushort_data_in;
            ushort_expected_fill_value = fv ? ushort_fill_value : NC_FILL_USHORT;
            fill_value = &ushort_expected_fill_value;
            fill_value_in = &ushort_fill_value_in;
            break;
         case NC_UINT:
            for (i = 0; i < DIMSIZE * DIMSIZE; i++)
               uint_data[i] = mpi_rank;
            data = uint_data;
            data_in = uint_data_in;
            uint_expected_fill_value = fv ? uint_fill_value : NC_FILL_UINT;
            fill_value = &uint_expected_fill_value;
            fill_value_in = &uint_fill_value_in;
            break;
         case NC_INT64:
            for (i = 0; i < DIMSIZE * DIMSIZE; i++)
               int64_data[i] = mpi_rank;
            data = int64_data;
            data_in = int64_data_in;
            int64_expected_fill_value = fv ? int64_fill_value : NC_FILL_INT64;
            fill_value = &int64_expected_fill_value;
            fill_value_in = &int64_fill_value_in;
            break;
         case NC_UINT64:
            for (i = 0; i < DIMSIZE * DIMSIZE; i++)
               uint64_data[i] = mpi_rank;
            data = uint64_data;
            data_in = uint64_data_in;
            uint64_expected_fill_value = fv ? uint64_fill_value : NC_FILL_UINT64;
            fill_value = &uint64_expected_fill_value;
            fill_value_in = &uint64_fill_value_in;
            break;
         }
      
         /* Create a file name. */
         sprintf(file_name, "%s_type_%d_fv_%d.nc", TEST_NAME, test_type[tt], fv);

         /* Create a parallel netcdf-4 file. */
         if (nc_create_par(file_name, NC_NETCDF4, comm, info, &ncid)) ERR;

         /* Get the type len. */
         if (nc_inq_type(ncid, test_type[tt], NULL, &type_size)) ERR;

         /* A global attribute holds the number of processors that created
          * the file. */
         if (nc_put_att_int(ncid, NC_GLOBAL, "num_processors", NC_INT, 1, &mpi_size)) ERR;

         /* Create three dimensions. */
         if (nc_def_dim(ncid, DIM1_NAME, NUM_SLABS, dimids)) ERR;
         if (nc_def_dim(ncid, DIM2_NAME, DIMSIZE, &dimids[1])) ERR;
         if (nc_def_dim(ncid, DIM3_NAME, DIMSIZE, &dimids[2])) ERR;

         /* Create one var. */
         if (nc_def_var(ncid, VAR_NAME, test_type[tt], NDIMS, dimids, &varid)) ERR;
         if (nc_put_att_int(ncid, varid, "var_num_processors", NC_INT, 1, &mpi_size)) ERR;
         if (fv == 1)
         {
            if (nc_def_var_fill(ncid, varid, NC_FILL, fill_value)) ERR;
            if (nc_inq_var_fill(ncid, varid, &fill_mode_in, fill_value_in)) ERR;
            if (fill_mode_in != NC_FILL) ERR;
            if (memcmp(fill_value_in, fill_value, type_size)) ERR;
         }
         else if (fv == 2)
         {
            if (nc_def_var_fill(ncid, varid, NC_NOFILL, NULL)) ERR;
            if (nc_inq_var_fill(ncid, varid, &fill_mode_in, NULL)) ERR;
            if (!fill_mode_in) ERR; /* nofill will be true */
         }            

         /* Write metadata to file. */
         if (nc_enddef(ncid)) ERR;

         /* Change access mode to collective, then back to independent. */
         if (nc_var_par_access(ncid, varid, NC_COLLECTIVE)) ERR;
         if (nc_var_par_access(ncid, varid, NC_INDEPENDENT)) ERR;

         if (!mpi_rank)
            start_time = MPI_Wtime();

         /* Write all the slabs this process is responsible for. */
         for (i = 0; i < NUM_SLABS / mpi_size; i++)
         {
            write_start[0] = NUM_SLABS / mpi_size * mpi_rank + i;

            /* Write one slab of data. Due to start/count settings,
             * every 16th value will be a fill value. */
            if (nc_put_vara(ncid, varid, write_start, write_count, data)) ERR;
         }

         /* On rank 0, keep track of time. */
         if (!mpi_rank)
         {
            total_time = MPI_Wtime() - start_time;
            printf("%d\t%g\t%g\n", mpi_size, total_time, DIMSIZE * DIMSIZE * NUM_SLABS *
                   sizeof(int) / total_time);
         }

         /* Close the netcdf file. */
         if (nc_close(ncid)) ERR;

         /* Reopen the file and check it. */
         if ((ret = nc_open_par(file_name, NC_NOWRITE, comm, info, &ncid))) ERR;
         if (nc_inq(ncid, &ndims_in, &nvars_in, &natts_in, &unlimdimid_in)) ERR;
         if (ndims_in != NDIMS || nvars_in != 1 || natts_in != 1 ||
             unlimdimid_in != -1) ERR;

         /* Check the attributes. */
         if (nc_get_att_int(ncid, NC_GLOBAL, "num_processors", &mpi_size_in)) ERR;
         if (mpi_size_in != mpi_size) ERR;
         if (nc_get_att_int(ncid, 0, "var_num_processors", &mpi_size_in)) ERR;
         if (mpi_size_in != mpi_size) ERR;
         if (fv == 1)
         {
            if (nc_inq_var_fill(ncid, varid, &fill_mode_in, fill_value_in)) ERR;
            if (fill_mode_in != NC_FILL) ERR;
            if (memcmp(fill_value_in, fill_value, type_size)) ERR;
         }

         /* Read all the slabs this process is responsible for. */
         for (i = 0; i < NUM_SLABS / mpi_size; i++)
         {
            read_start[0] = NUM_SLABS / mpi_size * mpi_rank + i;
            /* printf("mpi_rank %d i %d read_start[0] %ld\n", mpi_rank, i, read_start[0]); */

            /* Read one slab of data. */
            if (nc_get_vara(ncid, varid, read_start, read_count, data_in)) ERR;

            /* Check data.  For the third fill value test, fill is
             * turned off. So don't bother testing the values where k
             * is zero. */
            /* printf("mpi_rank %d fv %d i %d j %d k %d int_data_in[j * k] %d int_expected_fill_value %d " */
            /*        "expected_value %d\n", mpi_rank, fv, i, j, k, int_data_in[j * k], */
            /*        int_expected_fill_value, expected_value); */
            switch (test_type[tt])
            {
            case NC_BYTE:
               for (j = 0; j < DIMSIZE; j++)
                  for (k = 0; k < DIMSIZE; k++)
                     if (fv < 2 || k)
                        if (byte_data_in[j * DIMSIZE + k] != (signed char)(k ? mpi_rank : byte_expected_fill_value)) ERR;
               break;
            case NC_SHORT:
               for (j = 0; j < DIMSIZE; j++)
                  for (k = 0; k < DIMSIZE; k++)
                     if (fv < 2 || k)
                        if (short_data_in[j * DIMSIZE + k] != (short)(k ? mpi_rank : short_expected_fill_value)) ERR;
               break;
            case NC_INT:
               for (j = 0; j < DIMSIZE; j++)
                  for (k = 0; k < DIMSIZE; k++)
                     if (fv < 2 || k)
                        if (int_data_in[j * DIMSIZE + k] != (int)(k ? mpi_rank : int_expected_fill_value)) ERR;
               break;
            case NC_FLOAT:
               for (j = 0; j < DIMSIZE; j++)
                  for (k = 0; k < DIMSIZE; k++)
                     if (fv < 2 || k)
                        if (float_data_in[j * DIMSIZE + k] != (float)(k ? mpi_rank : float_expected_fill_value)) ERR;
               break;
            case NC_DOUBLE:
               for (j = 0; j < DIMSIZE; j++)
                  for (k = 0; k < DIMSIZE; k++)
                     if (fv < 2 || k)
                        if (double_data_in[j * DIMSIZE + k] != (double)(k ? mpi_rank : double_expected_fill_value)) ERR;
               break;
            case NC_UBYTE:
               for (j = 0; j < DIMSIZE; j++)
                  for (k = 0; k < DIMSIZE; k++)
                     if (fv < 2 || k)
                        if (ubyte_data_in[j * DIMSIZE + k] != (unsigned char)(k ? mpi_rank : ubyte_expected_fill_value)) ERR;
               break;
            case NC_USHORT:
               for (j = 0; j < DIMSIZE; j++)
                  for (k = 0; k < DIMSIZE; k++)
                     if (fv < 2 || k)
                        if (ushort_data_in[j * DIMSIZE + k] != (unsigned short)(k ? mpi_rank : ushort_expected_fill_value)) ERR;
               break;
            case NC_UINT:
               for (j = 0; j < DIMSIZE; j++)
                  for (k = 0; k < DIMSIZE; k++)
                     if (fv < 2 || k)
                        if (uint_data_in[j * DIMSIZE + k] != (unsigned int)(k ? mpi_rank : uint_expected_fill_value)) ERR;
               break;
            case NC_INT64:
               for (j = 0; j < DIMSIZE; j++)
                  for (k = 0; k < DIMSIZE; k++)
                     if (fv < 2 || k)
                        if (int64_data_in[j * DIMSIZE + k] != (long long int)(k ? mpi_rank : int64_expected_fill_value)) ERR;
               break;
            case NC_UINT64:
               for (j = 0; j < DIMSIZE; j++)
                  for (k = 0; k < DIMSIZE; k++)
                     if (fv < 2 || k)
                        if (uint64_data_in[j * DIMSIZE + k] != (unsigned long long int)(k ? mpi_rank : uint64_expected_fill_value)) ERR;
               break;
            }
         } /* next slab */

         /* Close the netcdf file. */
         if (nc_close(ncid))  ERR;

         if (!mpi_rank)
            SUMMARIZE_ERR;
      } /* next test type */
   } /* next fill value test run */
   
   /* Shut down MPI. */
   MPI_Finalize();

   if (!mpi_rank)
      FINAL_RESULTS;
   return 0;
}
