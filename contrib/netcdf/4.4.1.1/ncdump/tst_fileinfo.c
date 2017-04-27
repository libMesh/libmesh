/* This is part of the netCDF package. Copyright 2008 University
   Corporation for Atmospheric Research/Unidata See COPYRIGHT file for
   conditions of use. See www.unidata.ucar.edu for more info.
*/

/*
Test _NCProperties and other special attributes
*/

#include "config.h"
#ifdef HAVE_UNISTD_H
#include "unistd.h"
#endif

#include <hdf5.h>
#include "nc_tests.h"
#include "err_macros.h"
#include "netcdf.h"
#include "nc4internal.h"

#define NC4FILE "nc4_fileinfo.nc"
#define HDFFILE "hdf5_fileinfo.hdf"
#define INT_ATT_NAME "int_attr"
#define INT_VAR_NAME "int_var"
#define GROUPNAME "subgroup"
#define DIMNAME "time"

int
main(int argc, char **argv)
{
    printf("\n*** Testing 'Fileinfo attributes.\n");

    {
	hid_t fileid;
	hid_t fcplid;
	hid_t scalar_spaceid;

        printf("*** creating test file using HDF5 directly %s...", HDFFILE);

        /* Create scalar dataspace */
        if((scalar_spaceid = H5Screate(H5S_SCALAR)) < 0) ERR;

        /* Set creation ordering for file, so we can revise its contents later */
        if((fcplid = H5Pcreate(H5P_FILE_CREATE)) < 0) ERR;
        if(H5Pset_link_creation_order(fcplid, H5P_CRT_ORDER_TRACKED) < 0) ERR;
        if(H5Pset_attr_creation_order(fcplid, H5P_CRT_ORDER_TRACKED) < 0) ERR;

        /* Create new file, using default properties */
        if((fileid = H5Fcreate(HDFFILE, H5F_ACC_TRUNC, fcplid, H5P_DEFAULT)) < 0) ERR;
        /* Close file creation property list */
        if(H5Pclose(fcplid) < 0) ERR;

        /* Add attributes to root group */
        {
            hid_t scalar_spaceid = -1;
            hid_t attid = -1;

            /* Create scalar dataspace */
            if((scalar_spaceid = H5Screate(H5S_SCALAR)) < 0) ERR;

            /* Create attribute with native integer datatype on object */
            if((attid = H5Acreate2(fileid, INT_ATT_NAME, H5T_NATIVE_INT, scalar_spaceid, H5P_DEFAULT, H5P_DEFAULT)) < 0) ERR;
            if(H5Aclose(attid) < 0) ERR;

            /* Clean up objects created */
            if(H5Sclose(scalar_spaceid) < 0) ERR;
        }

        /* Close rest */
        if(H5Sclose(scalar_spaceid) < 0) ERR;
        if(H5Fclose(fileid) < 0) ERR;
    }

    {
	int root, grpid, varid, stat, natts, id;
	int data = 17;
	const char* sdata = "text";
	char ncprops[8192];
	size_t len;
	int dimid;
        nc_type xtype;
	char name[NC_MAX_NAME];

        printf("\n*** creating netcdf-4 test file using netCDF %s...", NC4FILE);

	if(nc_create(NC4FILE,NC_WRITE|NC_CLOBBER|NC_NETCDF4,&root)!=0) ERR;
	/* Create global attribute */
	if(nc_put_att_int(root,NC_GLOBAL,INT_ATT_NAME,NC_INT,1,&data)!=0) ERR;
	/* Create global variable */
	if(nc_def_var(root,INT_VAR_NAME,NC_INT,0,NULL,&varid)!=0) ERR;
	/* Create attribute on var */
	if(nc_put_att_int(root,varid,INT_ATT_NAME,NC_INT,1,&data)!=0) ERR;
	/* Create global subgroup */
	if(nc_def_grp(root,GROUPNAME,&grpid)!=0) ERR;
	/* Create global attribute in the group */
	if(nc_put_att_int(grpid,NC_GLOBAL,INT_ATT_NAME,NC_INT,1,&data)!=0) ERR;
	/* Create _NCProperties as var attr and as subgroup attribute */
	if(nc_put_att_text(grpid,NC_GLOBAL,NCPROPS,strlen(sdata),sdata)!=0) ERR;
	if(nc_put_att_text(root,varid,NCPROPS,strlen(sdata),sdata)!=0) ERR;
	/* Create var + dimension to cause e.g. dimscales to appear */
	if(nc_def_dim(root,DIMNAME,(size_t)4,&dimid)!=0) ERR;
	if(nc_def_var(root,DIMNAME,NC_INT,1,&dimid,&varid)!=0) ERR; /* same name */
	/* Close, then re-open */
	if(nc_close(root)) ERR;
	if(nc_open(NC4FILE,NC_WRITE|NC_NETCDF4,&root)!=0) ERR;

	/* Is all invisible attributes actually invisible vis-a-vis nc_inq? */
        if(nc_inq(root,NULL,NULL,&natts,NULL)!=0) ERR;
	if(natts != 1) ERR;

	/* Now, fiddle with the NCPROPS attribute */

	/* Get its metadata */
	if(nc_inq_att(root,NC_GLOBAL,NCPROPS,&xtype,&len)!=0) ERR;
	if(xtype != NC_CHAR) ERR;

	/* Read in two ways */
	if(nc_get_att_text(root,NC_GLOBAL,NCPROPS,ncprops)!=0) ERR;
	if(strlen(ncprops) != len) ERR;

	/* Attempt to get attribute metadata piecemeal; some will fail */
	id = -1;
	stat = nc_inq_attid(root,NC_GLOBAL,NCPROPS,&id);
	if(stat == NC_NOERR) ERR;
	stat = nc_inq_attname(root,NC_GLOBAL,id,name);
	if(stat == NC_NOERR) ERR;
	if(nc_inq_atttype(root,NC_GLOBAL,NCPROPS,&xtype)!=0) ERR;
	if(xtype != NC_CHAR) ERR;
	if(nc_inq_attlen(root,NC_GLOBAL,NCPROPS,&len)!=0) ERR;
	if(len != strlen(ncprops)) ERR;

	/*Overwrite _NCProperties root attribute; should fail */
	stat = nc_put_att_text(root,NC_GLOBAL,NCPROPS,strlen(sdata),sdata);
	if(stat == NC_NOERR) ERR;

	/* Delete; should fail */
	stat = nc_del_att(root,NC_GLOBAL,NCPROPS);
        if(stat != NC_ENOTATT) ERR;

	/* Ditto _SuperblockVersion */

	/* Get its metadata */
	if(nc_inq_att(root,NC_GLOBAL,SUPERBLOCKATT,&xtype,&len)!=0) ERR;
	if(xtype != NC_INT) ERR;
	if(len != 1) ERR;

	if(nc_get_att_int(root,NC_GLOBAL,SUPERBLOCKATT,&data)!=0) ERR;

	/* Attempt to get attribute metadata piecemeal */
	stat = nc_inq_attid(root,NC_GLOBAL,SUPERBLOCKATT,&id);
	if(stat == NC_NOERR) ERR;
	stat = nc_inq_attname(root,NC_GLOBAL,id,name);
	if(stat == NC_NOERR) ERR;
	if(nc_inq_atttype(root,NC_GLOBAL,SUPERBLOCKATT,&xtype)!=0) ERR;
	if(xtype != NC_INT) ERR;
	if(nc_inq_attlen(root,NC_GLOBAL,SUPERBLOCKATT,&len)!=0) ERR;
	if(len != 1) ERR;

	/*Overwrite; should fail */
	stat = nc_put_att_int(root,NC_GLOBAL,NCPROPS,NC_INT,1,&data);
	if(stat == NC_NOERR) ERR;

	/* Delete; should fail */
	stat = nc_del_att(root,NC_GLOBAL,SUPERBLOCKATT);
        if(stat == NC_NOERR) ERR;

	/* Ditto _IsNetcdf4 */

	/* Get its metadata */
	if(nc_inq_att(root,NC_GLOBAL,ISNETCDF4ATT,&xtype,&len)!=0) ERR;
	if(xtype != NC_INT) ERR;
	if(len != 1) ERR;

	if(nc_get_att_int(root,NC_GLOBAL,ISNETCDF4ATT,&data)!=0) ERR;

	/* Attempt to get attribute metadata piecemeal */
	stat = nc_inq_attid(root,NC_GLOBAL,ISNETCDF4ATT,&id);
	if(stat == NC_NOERR) ERR;
	stat = nc_inq_attname(root,NC_GLOBAL,id,name);
	if(stat == NC_NOERR) ERR;
	if(nc_inq_atttype(root,NC_GLOBAL,ISNETCDF4ATT,&xtype)!=0) ERR;
	if(xtype != NC_INT) ERR;
	if(nc_inq_attlen(root,NC_GLOBAL,ISNETCDF4ATT,&len)!=0) ERR;
	if(len != 1) ERR;

	/*Overwrite; should fail */
	stat = nc_put_att_int(root,NC_GLOBAL,ISNETCDF4ATT,NC_INT,1,&data);
	if(stat == NC_NOERR) ERR;

	/* Delete; should fail */
	stat = nc_del_att(root,NC_GLOBAL,ISNETCDF4ATT);
        if(stat == NC_NOERR) ERR;

	if(nc_close(root)!=0) ERR;
    }

    SUMMARIZE_ERR;
    FINAL_RESULTS;
}
