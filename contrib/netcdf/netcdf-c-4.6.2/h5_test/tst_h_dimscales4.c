/* This is part of the netCDF package.  Copyright 2005-2011 University
   Corporation for Atmospheric Research/Unidata See COPYRIGHT file for
   conditions of use.

   Test HDF5 file code. These are not intended to be exhaustive tests,
   but they use HDF5 the same way that netCDF-4 does, so if these
   tests don't work, than netCDF-4 won't work either.

*/

#include "h5_err_macros.h"
#include <hdf5.h>
#include <H5DSpublic.h>
#include <ncdimscale.h>

#define FILE_NAME "tst_h_dimscales4.h5"
#define DIMSCALE_NAME "AAA_dimscale"
#define VAR1_NAME "Watson"
#define VAR2_NAME "Holmes"
#define VAR3_NAME "Moriarty"
#define NDIMS 1
#define DIM_LEN 1
#define NAME_ATTRIBUTE "dimscale_name_attribute"
#define DIMSCALE_LABEL "dimscale_label"
#define STR_LEN 255
#define NC_MAX_NAME STR_LEN
#define NC_EHDFERR 255
#define DIM_WITHOUT_VARIABLE "This is a netCDF dimension but not a netCDF variable."

struct nc_hdf5_link_info 
{
   char name[STR_LEN];
   H5I_type_t obj_type;   
};   

herr_t alien_visitor(hid_t did, unsigned dim, hid_t dsid, 
		     void *visitor_data)
{
   H5G_stat_t statbuf;
   HDF5_OBJID_T *objid = visitor_data;

   /* Get more info on the dimscale object.*/
   if (H5Gget_objinfo(dsid, ".", 1, &statbuf) < 0) ERR;
   objid->fileno[0] = statbuf.fileno[0];
   objid->objno[0] = statbuf.objno[0];
   objid->fileno[1] = statbuf.fileno[1];
   objid->objno[1] = statbuf.objno[1];

   return 0;
}

int
main()
{
   printf("\n*** Checking HDF5 dimscales detach.\n");
   printf("*** Creating a file with two vars with one dimension scale...");
   {
      hid_t fileid, grpid, spaceid, var1_id, var2_id, dimscaleid;
      hid_t fcpl_id, fapl_id, create_propid, access_propid;
      hsize_t dims[NDIMS] = {DIM_LEN};
      char dimscale_wo_var[STR_LEN];
      float data = 42;

      /* Create file. */
      if ((fapl_id = H5Pcreate(H5P_FILE_ACCESS)) < 0) ERR;
      if (H5Pset_fclose_degree(fapl_id, H5F_CLOSE_STRONG)) ERR;
      if (H5Pset_cache(fapl_id, 0, CHUNK_CACHE_NELEMS, CHUNK_CACHE_SIZE,
		       CHUNK_CACHE_PREEMPTION) < 0) ERR;
      if (H5Pset_libver_bounds(fapl_id, H5F_LIBVER_LATEST, H5F_LIBVER_LATEST) < 0) ERR;
      if ((fcpl_id = H5Pcreate(H5P_FILE_CREATE)) < 0) ERR;
      if (H5Pset_link_creation_order(fcpl_id, (H5P_CRT_ORDER_TRACKED |
					       H5P_CRT_ORDER_INDEXED)) < 0) ERR;
      if (H5Pset_attr_creation_order(fcpl_id, (H5P_CRT_ORDER_TRACKED |
					       H5P_CRT_ORDER_INDEXED)) < 0) ERR;
      if ((fileid = H5Fcreate(FILE_NAME, H5F_ACC_TRUNC, fcpl_id, fapl_id)) < 0) ERR;
      if (H5Pclose(fapl_id) < 0) ERR;
      if (H5Pclose(fcpl_id) < 0) ERR;
      if ((grpid = H5Gopen2(fileid, "/", H5P_DEFAULT)) < 0) ERR;

      /* Create dimension scale. */
      if ((create_propid = H5Pcreate(H5P_DATASET_CREATE)) < 0) ERR;
      if (H5Pset_attr_creation_order(create_propid, H5P_CRT_ORDER_TRACKED|
				     H5P_CRT_ORDER_INDEXED) < 0) ERR;
      if ((spaceid = H5Screate_simple(1, dims, dims)) < 0) ERR;
      if ((dimscaleid = H5Dcreate1(grpid, DIMSCALE_NAME, H5T_IEEE_F32BE,
				   spaceid, create_propid)) < 0) ERR;
      if (H5Sclose(spaceid) < 0) ERR;
      if (H5Pclose(create_propid) < 0) ERR;
      sprintf(dimscale_wo_var, "%s%10d", DIM_WITHOUT_VARIABLE, DIM_LEN);
      if (H5DSset_scale(dimscaleid, dimscale_wo_var) < 0) ERR;

      /* Create a variable that uses this dimension scale. */
      if ((access_propid = H5Pcreate(H5P_DATASET_ACCESS)) < 0) ERR;
      if (H5Pset_chunk_cache(access_propid, CHUNK_CACHE_NELEMS,
   			     CHUNK_CACHE_SIZE, CHUNK_CACHE_PREEMPTION) < 0) ERR;
      if ((create_propid = H5Pcreate(H5P_DATASET_CREATE)) < 0) ERR;
      if (H5Pset_fill_value(create_propid, H5T_NATIVE_FLOAT, &data) < 0) ERR;
      if (H5Pset_layout(create_propid, H5D_CONTIGUOUS) < 0) ERR;
      if (H5Pset_attr_creation_order(create_propid, H5P_CRT_ORDER_TRACKED|
				     H5P_CRT_ORDER_INDEXED) < 0) ERR;
      if ((spaceid = H5Screate_simple(NDIMS, dims, dims)) < 0) ERR;
      if ((var1_id = H5Dcreate2(grpid, VAR1_NAME, H5T_NATIVE_FLOAT, spaceid, 
				H5P_DEFAULT, create_propid, access_propid)) < 0) ERR;
      if (H5Pclose(create_propid) < 0) ERR;
      if (H5Pclose(access_propid) < 0) ERR;
      if (H5Sclose(spaceid) < 0) ERR;
      if (H5DSattach_scale(var1_id, dimscaleid, 0) < 0) ERR;

      /* Create another variable that uses this dimension scale. */
      if ((access_propid = H5Pcreate(H5P_DATASET_ACCESS)) < 0) ERR;
      if (H5Pset_chunk_cache(access_propid, CHUNK_CACHE_NELEMS,
   			     CHUNK_CACHE_SIZE, CHUNK_CACHE_PREEMPTION) < 0) ERR;
      if ((create_propid = H5Pcreate(H5P_DATASET_CREATE)) < 0) ERR;
      if (H5Pset_fill_value(create_propid, H5T_NATIVE_FLOAT, &data) < 0) ERR;
      if (H5Pset_layout(create_propid, H5D_CONTIGUOUS) < 0) ERR;
      if (H5Pset_attr_creation_order(create_propid, H5P_CRT_ORDER_TRACKED|
				     H5P_CRT_ORDER_INDEXED) < 0) ERR;
      if ((spaceid = H5Screate_simple(NDIMS, dims, dims)) < 0) ERR;
      if ((var2_id = H5Dcreate2(grpid, VAR2_NAME, H5T_NATIVE_FLOAT, spaceid, 
				H5P_DEFAULT, create_propid, access_propid)) < 0) ERR;
      if (H5Pclose(create_propid) < 0) ERR;
      if (H5Pclose(access_propid) < 0) ERR;
      if (H5Sclose(spaceid) < 0) ERR;
      if (H5DSattach_scale(var2_id, dimscaleid, 0) < 0) ERR;

      /* Now detach the scales and remove the dimscale. This doesn't
       * work if I reverse the order of the statements. */
      if (H5DSdetach_scale(var2_id, dimscaleid, 0) < 0) ERR;     
      if (H5DSdetach_scale(var1_id, dimscaleid, 0) < 0) ERR;      

      /* Fold up our tents. */
      if (H5Dclose(var1_id) < 0) ERR;
      if (H5Dclose(dimscaleid) < 0) ERR;
      if (H5Gclose(grpid) < 0) ERR;
      if (H5Fclose(fileid) < 0) ERR;
   }
   SUMMARIZE_ERR;
   FINAL_RESULTS;
}

