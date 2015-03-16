/* This is part of the netCDF package.  Copyright 2013 University
   Corporation for Atmospheric Research/Unidata See COPYRIGHT file for
   conditions of use.

   This program sets up HDF5 files that contain scalar attributes and
   variables, of both string and numeric datatypes.  ncdump should handle
   all of these.
*/

#include <nc_tests.h>
#include <hdf5.h>

#define FILE_NAME "tst_h_scalar.nc"
#define VSTR_ATT1_NAME  "vstratt1"
#define VSTR_ATT2_NAME  "vstratt2"
#define VSTR_ATT3_NAME  "vstratt3"
#define VSTR_ATT4_NAME  "vstratt4"
#define VSTR_VAR1_NAME  "vstrvar1"
#define VSTR_VAR2_NAME  "vstrvar2"
#define VSTR_VAR3_NAME  "vstrvar3"
#define VSTR_VAR4_NAME  "vstrvar4"
#define FSTR_ATT_NAME   "fstratt"
#define FSTR_VAR_NAME   "fstrvar"
#define INT_ATT_NAME    "intatt"
#define INT_VAR_NAME    "intvar"

int
add_attrs(hid_t objid)
{
    hid_t scalar_spaceid = -1;
    hid_t vlstr_typeid = -1, fixstr_typeid = -1;
    char *vlstr;
    hid_t attid = -1;

    /* Create scalar dataspace */
    if ((scalar_spaceid = H5Screate(H5S_SCALAR)) < 0) ERR_GOTO;
    
    /* Create string datatypes */
    if ((vlstr_typeid = H5Tcreate(H5T_STRING, (size_t)H5T_VARIABLE)) < 0) ERR_GOTO;
    if ((fixstr_typeid = H5Tcreate(H5T_STRING, (size_t)10)) < 0) ERR_GOTO;


    /* Create attribute with VL string datatype on object */
    if ((attid = H5Acreate2(objid, VSTR_ATT1_NAME, vlstr_typeid, scalar_spaceid, H5P_DEFAULT, H5P_DEFAULT)) < 0) ERR_GOTO;
    /* No write, use fill value */
    if (H5Aclose(attid) < 0) ERR_GOTO;

    /* Create attribute with VL string datatype on object */
    if ((attid = H5Acreate2(objid, VSTR_ATT2_NAME, vlstr_typeid, scalar_spaceid, H5P_DEFAULT, H5P_DEFAULT)) < 0) ERR_GOTO;
    vlstr = NULL;
    if (H5Awrite(attid, vlstr_typeid, &vlstr) < 0) ERR_GOTO;
    if (H5Aclose(attid) < 0) ERR_GOTO;

    /* Create attribute with VL string datatype on object */
    if ((attid = H5Acreate2(objid, VSTR_ATT3_NAME, vlstr_typeid, scalar_spaceid, H5P_DEFAULT, H5P_DEFAULT)) < 0) ERR_GOTO;
    vlstr = malloc(10);
    *vlstr = '\0';
    if (H5Awrite(attid, vlstr_typeid, &vlstr) < 0) ERR_GOTO;
    if (H5Aclose(attid) < 0) ERR_GOTO;

    /* Create attribute with VL string datatype on object */
    if ((attid = H5Acreate2(objid, VSTR_ATT4_NAME, vlstr_typeid, scalar_spaceid, H5P_DEFAULT, H5P_DEFAULT)) < 0) ERR_GOTO;
    strcpy(vlstr, "foo");
    if (H5Awrite(attid, vlstr_typeid, &vlstr) < 0) ERR_GOTO;
    free(vlstr);
    if (H5Aclose(attid) < 0) ERR_GOTO;

    /* Create attribute with fixed-length string datatype on object */
    if ((attid = H5Acreate2(objid, FSTR_ATT_NAME, fixstr_typeid, scalar_spaceid, H5P_DEFAULT, H5P_DEFAULT)) < 0) ERR_GOTO;
    if (H5Aclose(attid) < 0) ERR_GOTO;

    /* Create attribute with native integer datatype on object */
    if ((attid = H5Acreate2(objid, INT_ATT_NAME, H5T_NATIVE_INT, scalar_spaceid, H5P_DEFAULT, H5P_DEFAULT)) < 0) ERR_GOTO;
    if (H5Aclose(attid) < 0) ERR_GOTO;


    /* Clean up objects created */
    if (H5Sclose(scalar_spaceid) < 0) ERR_GOTO;
    if (H5Tclose(vlstr_typeid) < 0) ERR_GOTO;
    if (H5Tclose(fixstr_typeid) < 0) ERR_GOTO;

    return(0);

error:
    H5E_BEGIN_TRY {
        H5Aclose(attid);
        H5Sclose(scalar_spaceid);
        H5Tclose(vlstr_typeid);
        H5Tclose(fixstr_typeid);
    } H5E_END_TRY;
    return(-1);
}

int
main()
{
    printf("\n*** Create file with datasets & attributes that have scalar dataspaces...");
    {
	hid_t fileid;
        hid_t fcplid;
	hid_t dsetid;
        hid_t dcplid;
	hid_t scalar_spaceid;
        hid_t vlstr_typeid, fixstr_typeid;
	hid_t attid;

        /* Create scalar dataspace */
	if ((scalar_spaceid = H5Screate(H5S_SCALAR)) < 0) ERR;

        /* Set creation ordering for file, so we can revise its contents later */
        if ((fcplid = H5Pcreate(H5P_FILE_CREATE)) < 0) ERR;
        if (H5Pset_link_creation_order(fcplid, H5P_CRT_ORDER_TRACKED) < 0) ERR;
        if (H5Pset_attr_creation_order(fcplid, H5P_CRT_ORDER_TRACKED) < 0) ERR;
	
	/* Create new file, using default properties */
	if ((fileid = H5Fcreate(FILE_NAME, H5F_ACC_TRUNC, fcplid, H5P_DEFAULT)) < 0) ERR;
	
        /* Close file creation property list */
        if (H5Pclose(fcplid) < 0) ERR;


        /* Create variable-length string datatype */
        if ((vlstr_typeid = H5Tcreate(H5T_STRING, (size_t)H5T_VARIABLE)) < 0) ERR;

        /* Create fixed-length string datatype */
        if ((fixstr_typeid = H5Tcreate(H5T_STRING, (size_t)10)) < 0) ERR;


        /* Set creation ordering for dataset, so we can revise its contents later */
        if ((dcplid = H5Pcreate(H5P_DATASET_CREATE)) < 0) ERR;
        if (H5Pset_attr_creation_order(dcplid, H5P_CRT_ORDER_TRACKED) < 0) ERR;

	
        /* Create scalar dataset with VL string datatype */
        if ((dsetid = H5Dcreate2(fileid, VSTR_VAR1_NAME, vlstr_typeid, scalar_spaceid, H5P_DEFAULT, dcplid, H5P_DEFAULT)) < 0) ERR;

        /* Add attributes to dataset */
        if (add_attrs(dsetid) < 0) ERR;

        /* Close VL string dataset */
        if (H5Dclose(dsetid) < 0) ERR;


        /* Create scalar dataset with fixed-length string datatype */
        if ((dsetid = H5Dcreate2(fileid, FSTR_VAR_NAME, fixstr_typeid, scalar_spaceid, H5P_DEFAULT, dcplid, H5P_DEFAULT)) < 0) ERR;

        /* Add attributes to dataset */
        if (add_attrs(dsetid) < 0) ERR;

        /* Close fixed-length string dataset */
        if (H5Dclose(dsetid) < 0) ERR;


        /* Create scalar dataset with native integer datatype */
        if ((dsetid = H5Dcreate2(fileid, INT_VAR_NAME, H5T_NATIVE_INT, scalar_spaceid, H5P_DEFAULT, dcplid, H5P_DEFAULT)) < 0) ERR;

        /* Add attributes to dataset */
        if (add_attrs(dsetid) < 0) ERR;

        /* Close native integer dataset */
        if (H5Dclose(dsetid) < 0) ERR;


        /* Add attributes to root group */
        if (add_attrs(fileid) < 0) ERR;


        /* Close dataset creation property list */
        if (H5Pclose(dcplid) < 0) ERR;

        /* Close string datatypes */
        if (H5Tclose(vlstr_typeid) < 0) ERR;
        if (H5Tclose(fixstr_typeid) < 0) ERR;


        /* Close rest */
	if (H5Sclose(scalar_spaceid) < 0) ERR;
	if (H5Fclose(fileid) < 0) ERR;
    }
    SUMMARIZE_ERR;

    printf("*** Revise file through netCDF-4 API...");
    {
	int ncid, varid;
        int ret;
        char *vlstr;

	if (nc_open(FILE_NAME, NC_WRITE, &ncid)) ERR;
        

        /* Define new VL string variable */
        if (nc_def_var(ncid, VSTR_VAR2_NAME , NC_STRING, 0, NULL, &varid)) ERR;

        /* Write to the variable */
        vlstr = NULL;
        if (nc_put_var(ncid, varid, &vlstr)) ERR;


        /* Define new VL string variable */
        if (nc_def_var(ncid, VSTR_VAR3_NAME , NC_STRING, 0, NULL, &varid)) ERR;

        /* Write to the variable */
        vlstr = malloc(10);
        *vlstr = '\0';
        if (nc_put_var(ncid, varid, &vlstr)) ERR;


        /* Define new VL string variable */
        if (nc_def_var(ncid, VSTR_VAR4_NAME , NC_STRING, 0, NULL, &varid)) ERR;

        /* Write to the variable */
        strcpy(vlstr, "foo");
        if (nc_put_var(ncid, varid, &vlstr)) ERR;
        free(vlstr);


	if (nc_close(ncid)) ERR;
    }
    SUMMARIZE_ERR;

    FINAL_RESULTS;
}

