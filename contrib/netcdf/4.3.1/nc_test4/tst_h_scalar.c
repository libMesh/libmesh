/* This is part of the netCDF package.  Copyright 2013 University
   Corporation for Atmospheric Research/Unidata See COPYRIGHT file for
   conditions of use.

   This program tests reading HDF5 files that contain scalar attributes and
   variables, of both string and numeric datatypes.  The netCDF-4 library
   should allow access to all of these.
*/

#include <nc_tests.h>
#include <hdf5.h>

#define FILE_NAME "tst_h_scalar.h5"
#define VSTR_ATT1_NAME  "vstratt1"
#define VSTR_ATT2_NAME  "vstratt2"
#define VSTR_ATT3_NAME  "vstratt3"
#define VSTR_ATT4_NAME  "vstratt4"
#define VSTR_VAR1_NAME  "vstrvar1"
#define VSTR_VAR2_NAME  "vstrvar2"
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
check_attrs(int ncid, int obj)
{
    int attid;
    int natts = 0;
    size_t len;
    nc_type type;
    char *vlstr;
    char fixstr[10];
    int x;

    /* Check the object's attributes are OK */
    if (nc_inq_varnatts(ncid, obj, &natts )) ERR_GOTO;
    if (natts != 6) ERR_GOTO;

    if (nc_inq_attid(ncid, obj, VSTR_ATT1_NAME, &attid)) ERR_GOTO;
    if (attid != 0) ERR_GOTO;
    if (nc_inq_atttype(ncid, obj, VSTR_ATT1_NAME, &type)) ERR_GOTO;
    if (type != NC_STRING) ERR_GOTO;
    if (nc_inq_attlen(ncid, obj, VSTR_ATT1_NAME, &len)) ERR_GOTO;
    if (len != 1) ERR_GOTO;
    vlstr = NULL;
    if (nc_get_att(ncid, obj, VSTR_ATT1_NAME, &vlstr)) ERR_GOTO;
    if (NULL != vlstr) ERR_GOTO;

    if (nc_inq_attid(ncid, obj, VSTR_ATT2_NAME, &attid)) ERR_GOTO;
    if (attid != 1) ERR_GOTO;
    if (nc_inq_atttype(ncid, obj, VSTR_ATT2_NAME, &type)) ERR_GOTO;
    if (type != NC_STRING) ERR_GOTO;
    if (nc_inq_attlen(ncid, obj, VSTR_ATT2_NAME, &len)) ERR_GOTO;
    if (len != 1) ERR_GOTO;
    vlstr = NULL;
    if (nc_get_att(ncid, obj, VSTR_ATT2_NAME, &vlstr)) ERR_GOTO;
    if (NULL != vlstr) ERR_GOTO;

    if (nc_inq_attid(ncid, obj, VSTR_ATT3_NAME, &attid)) ERR_GOTO;
    if (attid != 2) ERR_GOTO;
    if (nc_inq_atttype(ncid, obj, VSTR_ATT3_NAME, &type)) ERR_GOTO;
    if (type != NC_STRING) ERR_GOTO;
    if (nc_inq_attlen(ncid, obj, VSTR_ATT3_NAME, &len)) ERR_GOTO;
    if (len != 1) ERR_GOTO;
    vlstr = NULL;
    if (nc_get_att(ncid, obj, VSTR_ATT3_NAME, &vlstr)) ERR_GOTO;
    if (strcmp(vlstr, "")) ERR_GOTO;
    free(vlstr);

    if (nc_inq_attid(ncid, obj, VSTR_ATT4_NAME, &attid)) ERR_GOTO;
    if (attid != 3) ERR_GOTO;
    if (nc_inq_atttype(ncid, obj, VSTR_ATT4_NAME, &type)) ERR_GOTO;
    if (type != NC_STRING) ERR_GOTO;
    if (nc_inq_attlen(ncid, obj, VSTR_ATT4_NAME, &len)) ERR_GOTO;
    if (len != 1) ERR_GOTO;
    vlstr = NULL;
    if (nc_get_att(ncid, obj, VSTR_ATT4_NAME, &vlstr)) ERR_GOTO;
    if (strcmp(vlstr, "foo")) ERR_GOTO;
    free(vlstr);

    if (nc_inq_attid(ncid, obj, FSTR_ATT_NAME, &attid)) ERR_GOTO;
    if (attid != 4) ERR_GOTO;
    if (nc_inq_atttype(ncid, obj, FSTR_ATT_NAME, &type)) ERR_GOTO;
    if (type != NC_CHAR) ERR_GOTO;
    if (nc_inq_attlen(ncid, obj, FSTR_ATT_NAME, &len)) ERR_GOTO;
    if (len != 10) ERR_GOTO;
    memset(fixstr, 1, sizeof(fixstr));
    if (nc_get_att(ncid, obj, FSTR_ATT_NAME, fixstr)) ERR_GOTO;
    if ('\0' != fixstr[0]) ERR_GOTO;

    if (nc_inq_attid(ncid, obj, INT_ATT_NAME, &attid)) ERR_GOTO;
    if (attid != 5) ERR_GOTO;
    if (nc_inq_atttype(ncid, obj, INT_ATT_NAME, &type)) ERR_GOTO;
    if (type != NC_INT) ERR_GOTO;
    if (nc_inq_attlen(ncid, obj, INT_ATT_NAME, &len)) ERR_GOTO;
    if (len != 1) ERR_GOTO;
    x = -1;
    if (nc_get_att(ncid, obj, INT_ATT_NAME, &x)) ERR_GOTO;
    if (0 != x) ERR_GOTO;


    return(0);

error:
    return(-1);
}

int
main()
{
    printf("\n*** Creating file with datasets & attributes that have scalar dataspaces...");
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

    printf("*** Checking accessing file through netCDF-4 API...");
    {
	int ncid, varid;
        size_t len;
        nc_type type;
        int ndims;
        char *vlstr;
        int x;

	if (nc_open(FILE_NAME, NC_NOWRITE, &ncid)) ERR;

        /* Check the global attributes are OK */
        if (check_attrs(ncid, NC_GLOBAL) < 0) ERR;

        /* Verify that the VL string dataset is present and OK */
	if (nc_inq_varid(ncid, VSTR_VAR1_NAME, &varid)) ERR;
        if (varid != 0) ERR;
	if (nc_inq_vartype(ncid, varid, &type)) ERR;
	if (type != NC_STRING) ERR;
        if (nc_inq_varndims(ncid, varid, &ndims)) ERR;
        if (ndims != 0) ERR;
        vlstr = NULL;
        if (nc_get_var(ncid, varid, &vlstr)) ERR;
        if (NULL != vlstr) ERR;

        /* Check the variable's attributes are OK */
        if (check_attrs(ncid, varid) < 0) ERR;

        /* Verify that the fixed-length string dataset is present and OK */
	if (nc_inq_varid(ncid, FSTR_VAR_NAME, &varid)) ERR;
        if (varid != 1) ERR;
	if (nc_inq_vartype(ncid, varid, &type)) ERR;
	if (type != NC_STRING) ERR;
        if (nc_inq_varndims(ncid, varid, &ndims)) ERR;
        if (ndims != 0) ERR;
        vlstr = NULL;
        if (nc_get_var(ncid, varid, &vlstr)) ERR;
        if ('\0' != *vlstr) ERR;
        free(vlstr);

        /* Check the variable's attributes are OK */
        if (check_attrs(ncid, varid) < 0) ERR;

        /* Verify that the integer dataset is present and OK */
	if (nc_inq_varid(ncid, INT_VAR_NAME, &varid)) ERR;
        if (varid != 2) ERR;
	if (nc_inq_vartype(ncid, varid, &type)) ERR;
	if (type != NC_INT) ERR;
        if (nc_inq_varndims(ncid, varid, &ndims)) ERR;
        if (ndims != 0) ERR;
        x = -1;
        if (nc_get_var(ncid, varid, &x)) ERR;
        if (0 != x) ERR;

        /* Check the variable's attributes are OK */
        if (check_attrs(ncid, varid) < 0) ERR;

	if (nc_close(ncid)) ERR;
    }
    SUMMARIZE_ERR;

    printf("*** Checking revising file through netCDF-4 API...");
    {
	int ncid, varid;
        char *vlstr;

	if (nc_open(FILE_NAME, NC_WRITE, &ncid)) ERR;
        
        /* Write to the VL string variable */
	if (nc_inq_varid(ncid, VSTR_VAR1_NAME, &varid)) ERR;
        vlstr = NULL;
        if (nc_put_var(ncid, varid, &vlstr)) ERR;

        vlstr = malloc(10);
        *vlstr = '\0';
        if (nc_put_var(ncid, varid, &vlstr)) ERR;

        strcpy(vlstr, "foo");
        if (nc_put_var(ncid, varid, &vlstr)) ERR;
        free(vlstr);


        /* Write to a VL string attribute */
        vlstr = NULL;
        if (nc_put_att(ncid, varid, VSTR_ATT1_NAME, NC_STRING, 1, &vlstr)) ERR;

        vlstr = malloc(10);
        *vlstr = '\0';
        if (nc_put_att(ncid, varid, VSTR_ATT1_NAME, NC_STRING, 1, &vlstr)) ERR;

        strcpy(vlstr, "foo");
        if (nc_put_att(ncid, varid, VSTR_ATT1_NAME, NC_STRING, 1, &vlstr)) ERR;
        free(vlstr);


        /* Define a new VL string variable */
        if (nc_def_var(ncid, VSTR_VAR2_NAME , NC_STRING, 0, NULL, &varid)) ERR;

        /* Write to the variable's fill-value */
        vlstr = NULL;
        if (nc_put_att(ncid, varid, _FillValue, NC_STRING, 1, &vlstr)) ERR;

        vlstr = malloc(10);
        *vlstr = '\0';
        if (nc_put_att(ncid, varid, _FillValue, NC_STRING, 1, &vlstr)) ERR;

        strcpy(vlstr, "foo");
        if (nc_put_att(ncid, varid, _FillValue, NC_STRING, 1, &vlstr)) ERR;
        free(vlstr);


	if (nc_close(ncid)) ERR;
    }
    SUMMARIZE_ERR;

    FINAL_RESULTS;
}

