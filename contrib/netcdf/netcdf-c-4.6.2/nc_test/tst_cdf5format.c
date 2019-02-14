/* This is part of the netCDF package.
   Copyright 2005 University Corporation for Atmospheric Research/Unidata
   See COPYRIGHT file for conditions of use.

   Test fix of bug involving creation of a file with PnetCDF APIs,
   then opening and modifying the file with netcdf.

   Author: Wei-keng Liao, Ed Hartnett
*/
#include <config.h>
#include <nc_tests.h>
#include <err_macros.h>

#define NVARS 6
#define NX 5
#define NDIM2 2
#define FILENAME "tst_cdf5format.nc"

/* Write a file with 2 dims and 6 vars, including some sample data. */
int
write2(int ncid, int parallel)
{
   int dimid[NDIM2];
   char str[NC_MAX_NAME + 1];
   int varid[NVARS];
   
   /* define dimension */
   if (nc_def_dim(ncid, "Y", NC_UNLIMITED, &dimid[0])) ERR;
   if (nc_def_dim(ncid, "X", NX, &dimid[1])) ERR;

   /* Define vars. */
   for (int i = 0; i < NVARS; i++)
   {
      if (i % 2)
      {
	 sprintf(str, "fixed_var_%d",i);
	 if (nc_def_var(ncid, str, NC_INT, 1, &dimid[1], &varid[i])) ERR;
      }
      else
      {
	 sprintf(str, "record_var_%d",i);
	 if (nc_def_var(ncid, str, NC_INT, 2, dimid, &varid[i])) ERR;
      }
   }
   
   if (nc_enddef(ncid)) ERR;
   
   /* write all variables */
   for (int i = 0; i < NVARS; i++)
   {
      size_t start[NDIM2] = {0, 0};
      size_t count[NDIM2];
      int buf[NX];
      
      /* Initialize some data. */
      for (int j = 0; j < NX; j++)
	 buf[j] = i * 10 + j;

      /* Write the data. */
      if (i % 2)
      {
	 count[0] = NX;
	 if (nc_put_vara_int(ncid, varid[i], start, count, buf)) ERR;
      }
      else
      {
	 count[0] = 1;
	 count[1] = NX;
	 if (nc_put_vara_int(ncid, varid[i], start, count, buf)) ERR;
      }
   }
   return 0;
}

/* Add some attributes to the vars of an open file. */
int
extend(int ncid)
{
   int i;
   char str[32];

   if (nc_redef(ncid)) ERR;
   
   /* add attributes to make header grow */
   for (i = 0; i < NVARS; i++)
   {
      sprintf(str, "annotation_for_var_%d", i);
      if (nc_put_att_text(ncid, i, "text_attr", strlen(str), str)) ERR;
   }
   if (nc_enddef(ncid)) ERR;
   return NC_NOERR;
}

/* Read the file and check the data. */
int
read2(int ncid)
{
   for (int i = 0; i < NVARS; i++)
   {
      int buf[NX];
      size_t start[2] = {0, 0}, count[2];      

      if (i % 2)
      {
	 count[0] = NX;
      }
      else
      {
	 count[0] = 1;
	 count[1] = NX;
      }
      if (nc_get_vara_int(ncid, i, start, count, buf)) ERR;	 
      for (int j = 0; j < NX; j++)
      {
	 if (buf[j] != i * 10 + j)
	 {
	    printf("unexpected read value var i=%d buf[j=%d]=%d should be %d\n",
		   i, j, buf[j], i * 10 + j);
	    ERR;
	 }
      }
   }
   return 0;
}

int main(int argc, char* argv[])
{
   int rank, nprocs, ncid, cmode;
   MPI_Comm comm = MPI_COMM_SELF;
   MPI_Info info = MPI_INFO_NULL;

   MPI_Init(&argc, &argv);
   MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   if (rank > 0)
      return 2;

   printf("\nWrite using PNETCDF; Read using classic netCDF...");
   {
      /* Create a netCDF classic file with PnetCDF. */
      cmode = NC_CLOBBER;
      if (nc_create_par(FILENAME, cmode, comm, info, &ncid)) ERR;
      if (write2(ncid, 1)) ERR;
      if (nc_close(ncid)) ERR;
      
      /* Re-open the file with pnetCDF (parallel) and add var attributes. */
      if (nc_open_par(FILENAME, NC_WRITE, comm, info, &ncid)) ERR;
      if (extend(ncid)) ERR;
      if (nc_close(ncid)) ERR;
      
      /* Open with classic and check. */
      if (nc_open(FILENAME, 0, &ncid)) ERR;
      if (read2(ncid)) ERR;
      if (nc_close(ncid)) ERR;
   }
   SUMMARIZE_ERR;

   printf("\nWrite using CDF-5; Read using PNETCDF...");
   {
      /* Create a file with CDF5. */
      cmode = NC_CDF5 | NC_CLOBBER;
      if (nc_create(FILENAME, cmode, &ncid)) ERR;
      if (write2(ncid, 0)) ERR;
      if (nc_close(ncid)) ERR;
      
      /* Re-open the file with CDF5 and add some atts. */
      if (nc_open(FILENAME, NC_WRITE, &ncid)) ERR;
      if (extend(ncid)) ERR;
      if (nc_close(ncid)) ERR;

      /* Re-open with PnetCDF and check. */
      cmode = NC_NOCLOBBER;
      if (nc_open_par(FILENAME, cmode, comm, info, &ncid)) ERR;
      if (read2(ncid)) ERR;
      if (nc_close(ncid)) ERR;
   }
   SUMMARIZE_ERR;

   MPI_Finalize();
   FINAL_RESULTS;
   return 0;
}
