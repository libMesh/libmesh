/* This is part of the netCDF package.  Copyright 2005 University
   Corporation for Atmospheric Research/Unidata See COPYRIGHT file for
   conditions of use.

   This program tests fixes for reading netCDF-4 files that contain
   datasets with reference datatypes.  The netCDF-4 library should ignore
   the datasets & attributes that have reference datatypes and allow the 
   rest of the file to be accessed.
*/

#include <config.h>
#include <nc_tests.h>
#include <err_macros.h>
#include <hdf5.h>

#define FILE_NAME "tst_h_refs.h5"
#define REF_ATT_NAME "refatt"
#define REF_VAR_NAME "refvar"
#define INT_ATT_NAME "intatt"
#define INT_VAR_NAME "intvar"

int
main()
{
    printf("\n*** Creating file with datasets & attributes that have reference datatypes.\n");
    {
	hid_t fileid, scalar_spaceid;
	hid_t attid, dsetid;

	if ((scalar_spaceid = H5Screate(H5S_SCALAR)) < 0) ERR;
	
	/* Create new file, using default properties. */
	if ((fileid = H5Fcreate(FILE_NAME, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT)) < 0) ERR;
	

        /* Create dataset with reference datatype */
        if ((dsetid = H5Dcreate2(fileid, REF_VAR_NAME, H5T_STD_REF_OBJ, scalar_spaceid, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT)) < 0) ERR;

        /* Create attribute with reference datatype on reference dataset */
        if ((attid = H5Acreate2(dsetid, REF_ATT_NAME, H5T_STD_REF_OBJ, scalar_spaceid, H5P_DEFAULT, H5P_DEFAULT)) < 0) ERR;
        if (H5Aclose(attid) < 0) ERR;

        /* Create attribute with native int datatype on reference dataset */
        if ((attid = H5Acreate2(dsetid, INT_ATT_NAME, H5T_NATIVE_INT, scalar_spaceid, H5P_DEFAULT, H5P_DEFAULT)) < 0) ERR;
        if (H5Aclose(attid) < 0) ERR;

        /* Close reference dataset */
        if (H5Dclose(dsetid) < 0) ERR;


        /* Create dataset with native int datatype */
        if ((dsetid = H5Dcreate2(fileid, INT_VAR_NAME, H5T_NATIVE_INT, scalar_spaceid, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT)) < 0) ERR;

        /* Create attribute with reference datatype on integer dataset */
        if ((attid = H5Acreate2(dsetid, REF_ATT_NAME, H5T_STD_REF_OBJ, scalar_spaceid, H5P_DEFAULT, H5P_DEFAULT)) < 0) ERR;
        if (H5Aclose(attid) < 0) ERR;

        /* Create attribute with native int datatype on integer dataset */
        if ((attid = H5Acreate2(dsetid, INT_ATT_NAME, H5T_NATIVE_INT, scalar_spaceid, H5P_DEFAULT, H5P_DEFAULT)) < 0) ERR;
        if (H5Aclose(attid) < 0) ERR;

        /* Close integer dataset */
        if (H5Dclose(dsetid) < 0) ERR;


        /* Create attribute on root group with reference datatype */
        if ((attid = H5Acreate2(fileid, REF_ATT_NAME, H5T_STD_REF_OBJ, scalar_spaceid, H5P_DEFAULT, H5P_DEFAULT)) < 0) ERR;
        if (H5Aclose(attid) < 0) ERR;

        /* Create attribute on root group with native int datatype */
        if ((attid = H5Acreate2(fileid, INT_ATT_NAME, H5T_NATIVE_INT, scalar_spaceid, H5P_DEFAULT, H5P_DEFAULT)) < 0) ERR;
        if (H5Aclose(attid) < 0) ERR;


        /* Close rest */
	if (H5Sclose(scalar_spaceid) < 0) ERR;
	if (H5Fclose(fileid) < 0) ERR;
    }

    printf("*** Checking accessing file through netCDF-4 API...");
    {
	int ncid, varid, attid;
        int natts = 0;
	nc_type type;

	if (nc_open(FILE_NAME, NC_NOWRITE, &ncid)) ERR;

        /* Check the root group's attributes are OK */
	if (nc_inq_varnatts(ncid, NC_GLOBAL, &natts )) ERR;
	if (natts != 1) ERR;    /* Reference attribute should not be present */
	if (nc_inq_attid(ncid, NC_GLOBAL, INT_ATT_NAME, &attid)) ERR;
	if (nc_inq_atttype(ncid, NC_GLOBAL, INT_ATT_NAME, &type)) ERR;
	if (type != NC_INT) ERR;

        /* Verify that the reference dataset is not present */
	if (!nc_inq_varid(ncid, REF_VAR_NAME, &varid)) ERR;

        /* Verify that the integer dataset is present and OK */
	if (nc_inq_varid(ncid, INT_VAR_NAME, &varid)) ERR;
	if (nc_inq_vartype(ncid, varid, &type)) ERR;
	if (type != NC_INT) ERR;

        /* Check the integer dataset's attributes are OK */
	if (nc_inq_varnatts(ncid, varid, &natts )) ERR;
	if (natts != 1) ERR;    /* Reference attribute should not be present */
	if (nc_inq_attid(ncid, varid, INT_ATT_NAME, &attid)) ERR;
	if (nc_inq_atttype(ncid, varid, INT_ATT_NAME, &type)) ERR;
	if (type != NC_INT) ERR;

	if (nc_close(ncid)) ERR;
    }
    SUMMARIZE_ERR;

    FINAL_RESULTS;
}

