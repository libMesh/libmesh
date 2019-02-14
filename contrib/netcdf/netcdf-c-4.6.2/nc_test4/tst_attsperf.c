/* This is part of the netCDF package. Copyright 2018 University
 * Corporation for Atmospheric Research/Unidata. See COPYRIGHT file
 * for conditions of use.
 *
 * Test the netCDF-4 attribute code.
 *
 * WARNING: do not attempt to run this under windows because of the use
 * of gettimeofday().
 *
 * Ed Hartnett 6/19/18
*/

#include <config.h>
#include <nc_tests.h>
#include "err_macros.h"
#include "nc4internal.h"
#include <sys/time.h>
#include <hdf5.h>

#define TEST "tst_attsperf"
#define VAR "bigvar"
#define NDIMS 2
#define DIM0 "d0"
#define DIM1 "d1"
#define DIMSIZE0 16
#define DIMSIZE1 512
#define TOTALSIZE (DIMSIZE0 * DIMSIZE1)
#define NUM_ATTS 100
#define ATT_LEN 100
#define NUM_VARS 1

int
add_attributes(int ncid, int varid, size_t num_atts, size_t att_len)
{
   char att_name[NC_MAX_NAME + 1];
   double *att_data;
   int i, a;

   /* Allocate space for attribute data. */
   if (!(att_data = malloc(att_len * sizeof(double))))
      return NC_ENOMEM;

   /* Fill up data. */
   for (i = 0; i < ATT_LEN; i++)
      att_data[i] = i;

   /* Write a bunch of attributes. */
   for (a = 0; a < num_atts; a++)
   {
      sprintf(att_name, "%s_varid_%d_att_%d", TEST, varid, a);
      if (nc_put_att_double(ncid, varid, att_name, NC_DOUBLE,
                            att_len, att_data)) ERR;
   }

   free(att_data);

   return 0;
}

/* Build the test file. */
int
buildfile(size_t num_vars, size_t num_atts, size_t att_len,
          char *file_name)
{
   int ncid, varid;
   int dimids[NDIMS];
   int v;

   if (nc_create(file_name, NC_NETCDF4, &ncid)) ERR;

   if (nc_def_dim(ncid, DIM0, DIMSIZE0, &dimids[0])) ERR;
   if (nc_def_dim(ncid, DIM1, DIMSIZE1, &dimids[1])) ERR;
   for (v = 0; v < num_vars; v++)
   {
      char var_name[NC_MAX_NAME + 1];
      sprintf(var_name, "%s_var_%d", TEST, v);
      if (nc_def_var(ncid, var_name, NC_INT, NDIMS, dimids, &varid)) ERR;
      if (add_attributes(ncid, v, num_atts, att_len)) ERR;
   }
   if (!num_vars)
      if (add_attributes(ncid, NC_GLOBAL, num_atts, att_len)) ERR;
   if (nc_enddef(ncid)) ERR;

   if (nc_close(ncid)) ERR;
   return 0;
}

/* Open/close the file with netCDF. */
int
readfile(char *file_name, long long *delta, int do_inq, int num_vars)
{
   int ncid;
   struct timeval starttime, endtime;
   long long startt, endt;
   int natts;
   int v;

   /* Start the clock. */
   gettimeofday(&starttime, NULL);

   /* Open the file. */
   if (nc_open(file_name, NC_NETCDF4, &ncid)) ERR;

   /* Do an inq if desired, triggering read of atts. */
   for (v = 0; v < num_vars; v++)
      if (nc_inq_varnatts(ncid, v, &natts)) ERR;

   if (nc_inq_natts(ncid, &natts)) ERR;

   /* Close the file. */
   if (nc_close(ncid)) ERR;

   gettimeofday(&endtime, NULL);

   /* Compute the time delta */
   startt = (1000000 * starttime.tv_sec) + starttime.tv_usec;
   endt = (1000000 * endtime.tv_sec) + endtime.tv_usec;
   *delta = endt - startt;

   return 0;
}

/* Open/close the file with HDF5. */
int
readfile_hdf5(char *file_name, long long *delta, int do_inq, int num_vars)
{
   hid_t hdfid, hdf_grpid;
   hid_t fapl_id;
   struct timeval starttime, endtime;
   long long startt, endt;

   /* Start the clock. */
   gettimeofday(&starttime, NULL);

   /* Open and close the root group. */
   if ((fapl_id = H5Pcreate(H5P_FILE_ACCESS)) < 0) ERR;
   if (H5Pset_fclose_degree(fapl_id, H5F_CLOSE_SEMI)) ERR;
   if ((hdfid = H5Fopen(file_name, H5F_ACC_RDONLY, fapl_id)) < 0) ERR;
   if ((hdf_grpid = H5Gopen2(hdfid, "/", H5P_DEFAULT)) < 0) ERR;

   /* Do we want to do an inq? */
   if (do_inq)
   {
      if (num_vars)
      {
      }
      else /* global atts */
      {
         hsize_t num_obj;

         /* Find out how many attributes. */
         if ((num_obj = H5Aget_num_attrs(hdf_grpid)) < 0) ERR;
      }
   }

   if (H5Gclose(hdf_grpid) < 0) ERR;
   if (H5Fclose(hdfid) < 0) ERR;

   gettimeofday(&endtime, NULL);

   /* Compute the time delta */
   startt = (1000000 * starttime.tv_sec) + starttime.tv_usec;
   endt = (1000000 * endtime.tv_sec) + endtime.tv_usec;
   *delta = endt - startt;

   return 0;
}

#define NUM_RUNS 5
#define NUM_STEPS 20
#define FACTOR 100
#define NUM_INQ_TESTS 2
int
main(int argc, char **argv)
{
   size_t num_atts = 1;
   char file_name[NC_MAX_NAME + 1];
   float tot_nc4, tot_hdf5;
   int factor;
   int r, s, num_vars, do_inq;

   for (do_inq = 0; do_inq < NUM_INQ_TESTS; do_inq++)
   {
      for (num_vars = 0; num_vars <= NUM_VARS; num_vars++)
      {
         /* Reset. */
         num_atts = 1;

         factor = FACTOR;

         printf("*** %s %s\n", num_vars ? "variable attributes" : "global attributes",
                do_inq ? "with inq" : "");
         printf("Number of Attributes\tHDF5 Open Time (s)\tNetcdf4 Open Time (s)\n");
         for (s = 0; s < NUM_STEPS; s++)
         {
            tot_nc4 = 0;
            tot_hdf5 = 0;
            num_atts += factor * s;

            for (r = 0; r < NUM_RUNS; r++)
            {
               long long nc4_open_time;
               long long hdf5_open_time;

               /* Determine file name. */
               sprintf(file_name, "%s_%d_%d_%d.nc", TEST, num_vars, s, r);

               if (buildfile(num_vars, num_atts, ATT_LEN, file_name)) ERR;
               if (readfile(file_name, &nc4_open_time, do_inq, num_vars)) ERR;
               if (readfile_hdf5(file_name, &hdf5_open_time, do_inq, num_vars)) ERR;
               tot_nc4 += nc4_open_time;
               tot_hdf5 += hdf5_open_time;
            }

            /* Print average results to the millisec */
            printf("%ld\t%g\t%g\n", num_atts, tot_hdf5/((float)NUM_RUNS * 1000000),
                   tot_nc4/((float)NUM_RUNS * 1000000));
         }
      }
   } /* next do_inq */
   SUMMARIZE_ERR;
   FINAL_RESULTS;
}
