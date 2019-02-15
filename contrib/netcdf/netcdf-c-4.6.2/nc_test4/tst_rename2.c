/*
 * Test more renames of vars and dims.
 *
 * Ed Hartnett
 */

#include "nc_tests.h"
#include "err_macros.h"

#define TEST_NAME "tst_rename"
#define LAT "lat"
#define LON "lon"
#define LEV "lev"
#define DIM_X "x"
#define DIM_Y "y"
#define DIM_Z "z"

#define DIM1_LEN 4
#define NDIM1 1
#define NDIM3 3
#define NUM_ENDDEF_SETTINGS 2

int
main(int argc, char **argv)
{
#define NUM_FORMATS 2
   int formats[NUM_FORMATS] = {NC_FORMAT_NETCDF4, NC_FORMAT_NETCDF4_CLASSIC};
   int format;

   fprintf(stderr,"*** Testing more renames\n");

   for (format = 0; format < NUM_FORMATS; format++)
   {
      fprintf(stderr,"*** test renaming 3 dimensions with format %d...",
              formats[format]);
      {
         char filename[NC_MAX_NAME + 1];
         int ncid, dimid[NDIM3];
         int dimid_in;
         int enddef_setting;

         if (nc_set_default_format(formats[format], NULL)) ERR;

         for (enddef_setting = 0; enddef_setting < NUM_ENDDEF_SETTINGS;
              enddef_setting++)
         {
            sprintf(filename, "%s_%d_%d.nc", TEST_NAME, formats[format],
                    enddef_setting);
            
            /* Create file with three dims. */
            if (nc_create(filename, 0, &ncid)) ERR;
            if (nc_def_dim(ncid, LAT, DIM1_LEN, &dimid[0])) ERR;
            if (nc_def_dim(ncid, LON, DIM1_LEN, &dimid[1])) ERR;
            if (nc_def_dim(ncid, LEV, DIM1_LEN, &dimid[2])) ERR;

            if (enddef_setting)
            {
               if (nc_enddef(ncid)) ERR;
               if (nc_redef(ncid)) ERR;
            }
            
            /* Rename the dimensions. */
            if (nc_rename_dim(ncid, 0, DIM_X)) ERR;
            if (nc_rename_dim(ncid, 1, DIM_Y)) ERR;
            if (nc_rename_dim(ncid, 2, DIM_Z)) ERR;
            
            /* Close the file. */
            if (nc_close(ncid)) ERR;
            
            /* Reopen the file and check. */
            if (nc_open(filename, NC_NOWRITE, &ncid)) ERR;
            if (nc_inq_dimid(ncid, DIM_X, &dimid_in)) ERR;
            if (dimid_in != 0) ERR;
            if (nc_inq_dimid(ncid, DIM_Y, &dimid_in)) ERR;
            if (dimid_in != 1) ERR;
            if (nc_inq_dimid(ncid, DIM_Z, &dimid_in)) ERR;
            if (dimid_in != 2) ERR;
            if (nc_close(ncid)) ERR;
         } /* next enddef setting */
      }
      SUMMARIZE_ERR;

      fprintf(stderr,"*** test renaming 3 dims with coord data format %d...",
              formats[format]);
      {
         char filename[NC_MAX_NAME + 1];
         int ncid, dimid[NDIM3], varid[NDIM3];
         int dimid_in;
         int lat_data[DIM1_LEN] = {0, 1, 2, 3};
         int lon_data[DIM1_LEN] = {0, 10, 20, 30};
         int lev_data[DIM1_LEN] = {0, 100, 200, 300};

         if (nc_set_default_format(formats[format], NULL)) ERR;

         sprintf(filename, "%s_data_%d.nc", TEST_NAME, formats[format]);
         
         /* Create file with three dims. */
         if (nc_create(filename, 0, &ncid)) ERR;
         if (nc_def_dim(ncid, LAT, DIM1_LEN, &dimid[0])) ERR;
         if (nc_def_dim(ncid, LON, DIM1_LEN, &dimid[1])) ERR;
         if (nc_def_dim(ncid, LEV, DIM1_LEN, &dimid[2])) ERR;

         /* Define coordinate data vars. */
         if (nc_def_var(ncid, LAT, NC_INT, NDIM1, &dimid[0], &varid[0])) ERR;
         if (nc_def_var(ncid, LON, NC_INT, NDIM1, &dimid[1], &varid[1])) ERR;
         if (nc_def_var(ncid, LEV, NC_INT, NDIM1, &dimid[2], &varid[2])) ERR;

         if (nc_enddef(ncid)) ERR;

         if (nc_put_var(ncid, 0, lat_data)) ERR;
         if (nc_put_var(ncid, 1, lon_data)) ERR;
         if (nc_put_var(ncid, 2, lev_data)) ERR;

         if (nc_close(ncid)) ERR;
         if (nc_open(filename, NC_WRITE, &ncid)) ERR;
         if (nc_redef(ncid)) ERR;
         
         /* Rename the dimensions. */
         if (nc_rename_dim(ncid, 0, DIM_X)) ERR;
         if (nc_rename_dim(ncid, 1, DIM_Y)) ERR;
         if (nc_rename_dim(ncid, 2, DIM_Z)) ERR;
         
         /* Close the file. */
         if (nc_close(ncid)) ERR;
         
         /* Reopen the file and check. */
         if (nc_open(filename, NC_NOWRITE, &ncid)) ERR;
         if (nc_inq_dimid(ncid, DIM_X, &dimid_in)) ERR;
         if (dimid_in != 0) ERR;
         if (nc_inq_dimid(ncid, DIM_Y, &dimid_in)) ERR;
         if (dimid_in != 1) ERR;
         if (nc_inq_dimid(ncid, DIM_Z, &dimid_in)) ERR;
         if (dimid_in != 2) ERR;
         if (nc_close(ncid)) ERR;
      }
      SUMMARIZE_ERR;
      
   } /* next format */
   FINAL_RESULTS;
}
