/** \file 
Dimension functions

These functions define and inquire about dimensions.

Copyright 2010 University Corporation for Atmospheric
Research/Unidata. See COPYRIGHT file for more info.
*/

#include "ncdispatch.h"

/*! \defgroup dimensions Dimensions 

Dimensions are used to define the shape of data in netCDF.

Dimensions for a netCDF dataset are defined when it is created, while
the netCDF dataset is in define mode. Additional dimensions may be
added later by reentering define mode. A netCDF dimension has a name
and a length. In a netCDF classic or 64-bit offset file, at most one
dimension can have the unlimited length, which means variables using
this dimension can grow along this dimension. In a netCDF-4 file
multiple unlimited dimensions are supported.

There is a suggested limit (100) to the number of dimensions that can
be defined in a single netCDF dataset. The limit is the value of the
predefined macro NC_MAX_DIMS. The purpose of the limit is to make
writing generic applications simpler. They need only provide an array
of NC_MAX_DIMS dimensions to handle any netCDF dataset. The
implementation of the netCDF library does not enforce this advisory
maximum, so it is possible to use more dimensions, if necessary, but
netCDF utilities that assume the advisory maximums may not be able to
handle the resulting netCDF datasets.

Ordinarily, the name and length of a dimension are fixed when the
dimension is first defined. The name may be changed later, but the
length of a dimension (other than the unlimited dimension) cannot be
changed without copying all the data to a new netCDF dataset with a
redefined dimension length.

Dimension lengths in the C interface are type size_t rather than type
int to make it possible to access all the data in a netCDF dataset on
a platform that only supports a 16-bit int data type, for example
MSDOS. If dimension lengths were type int instead, it would not be
possible to access data from variables with a dimension length greater
than a 16-bit int can accommodate.

A netCDF dimension in an open netCDF dataset is referred to by a small
integer called a dimension ID. In the C interface, dimension IDs are
0, 1, 2, ..., in the order in which the dimensions were defined.

Operations supported on dimensions are:
- Create a dimension, given its name and length.
- Get a dimension ID from its name.
- Get a dimension's name and length from its ID.
- Rename a dimension. 
*/
/**@{*/

/*!  Define a new dimension. The function nc_def_dim adds a new
dimension to an open netCDF dataset in define mode. It returns (as an
argument) a dimension ID, given the netCDF ID, the dimension name, and
the dimension length. At most one unlimited length dimension, called
the record dimension, may be defined for each classic or 64-bit offset
netCDF dataset. NetCDF-4 datasets may have multiple unlimited
dimensions.

\param ncid NetCDF or group ID, from a previous call to nc_open(),
nc_create(), nc_def_grp(), or associated inquiry functions such as 
nc_inq_ncid().

\param name Name of the dimension to be created.

\param len Length of the dimension to be created. Use NC_UNLIMITED for
unlimited dimensions.

\param idp Pointer where dimension ID will be stored.

\retval ::NC_NOERR No error.
\returns ::NC_EBADID Not a valid ID.
\returns ::NC_ENOTINDEFINE Not in define mode.
\returns ::NC_EDIMSIZE Invalid dimension size.
\returns ::NC_EUNLIMIT NC_UNLIMITED size already in use
\returns ::NC_EMAXDIMS NC_MAX_DIMS exceeded
\returns ::NC_ENAMEINUSE String match to name in use
\returns ::NC_ENOMEM Memory allocation (malloc) failure
\returns ::NC_EPERM Write to read only

\section Example

Here is an example using nc_def_dim() to create a dimension named lat of
length 18 and a unlimited dimension named rec in a new netCDF dataset
named foo.nc:

\code
     #include <netcdf.h>
        ...
     int status, ncid, latid, recid;
        ...
     status = nc_create("foo.nc", NC_NOCLOBBER, &ncid);
     if (status != NC_NOERR) handle_error(status);
        ...
     status = nc_def_dim(ncid, "lat", 18L, &latid);
     if (status != NC_NOERR) handle_error(status);
     status = nc_def_dim(ncid, "rec", NC_UNLIMITED, &recid);
     if (status != NC_NOERR) handle_error(status);
\endcode

 */
int
nc_def_dim(int ncid, const char *name, size_t len, int *idp)
{
    NC* ncp;
    int stat = NC_check_id(ncid, &ncp);
    if(stat != NC_NOERR) return stat;
    return ncp->dispatch->def_dim(ncid, name, len, idp);
}

/*!
Find the ID of a dimension from the name.

The function nc_inq_dimid returns (as an argument) the ID of a netCDF
dimension, given the name of the dimension. If ndims is the number of
dimensions defined for a netCDF dataset, each dimension has an ID
between 0 and ndims-1.

\param ncid NetCDF or group ID, from a previous call to nc_open(),
nc_create(), nc_def_grp(), or associated inquiry functions such as 
nc_inq_ncid().

\param name Name of the dimension.

\param idp Pointer where dimension ID will be stored.

\returns ::NC_NOERR   No error.
\returns ::NC_EBADID  Not a valid ID.
\returns ::NC_EBADDIM Invalid dimension ID or name.
 */
int
nc_inq_dimid(int ncid, const char *name, int *idp)
{
    NC* ncp;
    int stat = NC_check_id(ncid, &ncp);
    if(stat != NC_NOERR) return stat;
    return ncp->dispatch->inq_dimid(ncid,name,idp);
}

/*!
Find the name and length of a dimension.

The length for the unlimited dimension, if any, is the number of
records written so far.

\param ncid NetCDF or group ID, from a previous call to nc_open(),
nc_create(), nc_def_grp(), or associated inquiry functions such as 
nc_inq_ncid().

\param dimid Dimension ID, from a previous call to nc_inq_dimid() or
nc_def_dim().

\param name Returned dimension name. The caller must allocate space
for the returned name. The maximum possible length, in characters, of
a dimension name is given by the predefined constant
NC_MAX_NAME. (This doesn't include the null terminator, so declare
your array to be size NC_MAX_NAME+1). The returned character array
will be null-terminated.

\param lenp Pointer to location for returned length of dimension. For
the unlimited dimension, this is the number of records written so far.

\returns ::NC_NOERR   No error.
\returns ::NC_EBADID  Not a valid ID.
\returns ::NC_EBADDIM Invalid dimension ID or name.

\section Example

Here is an example using nc_inq_dim() to determine the length of a
dimension named lat, and the name and current maximum length of the
unlimited dimension for an existing netCDF dataset named foo.nc:

\code
     #include <netcdf.h>
        ...
     int status, ncid, latid, recid;
     size_t latlength, recs;
     char recname[NC_MAX_NAME+1];
        ...
     status = nc_open("foo.nc", NC_NOWRITE, &ncid); 
     if (status != NC_NOERR) handle_error(status);
     status = nc_inq_unlimdim(ncid, &recid);
     if (status != NC_NOERR) handle_error(status);
        ...
     status = nc_inq_dimid(ncid, "lat", &latid);
     if (status != NC_NOERR) handle_error(status);
     status = nc_inq_dimlen(ncid, latid, &latlength);
     if (status != NC_NOERR) handle_error(status);

     status = nc_inq_dim(ncid, recid, recname, &recs);
     if (status != NC_NOERR) handle_error(status);
\endcode
 */
int
nc_inq_dim(int ncid, int dimid, char *name, size_t *lenp)
{
    NC* ncp;
    int stat = NC_check_id(ncid, &ncp);
    if(stat != NC_NOERR) return stat;
    return ncp->dispatch->inq_dim(ncid,dimid,name,lenp);
}

/*!
Rename a dimension.

This function renames an existing dimension in a netCDF dataset open
for writing. You cannot rename a dimension to have the same name as
another dimension.

For netCDF classic and 64-bit offset files, if the new name is longer
than the old name, the netCDF dataset must be in define mode.

For netCDF-4 files the dataset is switched to define more for the
rename, regardless of the name length.

\param ncid NetCDF or group ID, from a previous call to nc_open(),
nc_create(), nc_def_grp(), or associated inquiry functions such as 
nc_inq_ncid().

\param dimid Dimension ID, from a previous call to nc_inq_dimid() or
nc_def_dim().

\param name New name for dimension. Must be a null-terminated string
with length less than NC_MAX_NAME.

\returns ::NC_NOERR      No error.
\returns ::NC_EBADID     Not a valid ID.
\returns ::NC_EBADDIM    Invalid dimension ID or name.
\returns ::NC_ENAMEINUSE String match to name in use
\returns ::NC_ENOMEM     Memory allocation (malloc) failure
\returns ::NC_EPERM      Write to read only
\section Example

Here is an example using nc_rename_dim to rename the dimension lat to
latitude in an existing netCDF dataset named foo.nc:

\code
     #include <netcdf.h>
        ...
     int status, ncid, latid;
        ...
     status = nc_open("foo.nc", NC_WRITE, &ncid); 
     if (status != NC_NOERR) handle_error(status);
        ...
     status = nc_redef(ncid); 
     if (status != NC_NOERR) handle_error(status);
     status = nc_inq_dimid(ncid, "lat", &latid);
     if (status != NC_NOERR) handle_error(status);
     status = nc_rename_dim(ncid, latid, "latitude");
     if (status != NC_NOERR) handle_error(status);
     status = nc_enddef(ncid);
     if (status != NC_NOERR) handle_error(status);
\endcode
 */
int
nc_rename_dim(int ncid, int dimid, const char *name)
{
    NC* ncp;
    int stat = NC_check_id(ncid, &ncp);
    if(stat != NC_NOERR) return stat;
    return ncp->dispatch->rename_dim(ncid,dimid,name);
}

/*!
Find the number of dimensions.

In a classic model netCDF file, this funtion returns the number of
defined dimensions. In a netCDF-4/HDF5 file, this function returns the
number of dimensions available in the group specified by ncid, which
may be less than the total number of dimensions in a file. In a
netCDF-4/HDF5 file, dimensions are in all sub-groups, sub-sub-groups,
etc.

\param ncid NetCDF or group ID, from a previous call to nc_open(),
nc_create(), nc_def_grp(), or associated inquiry functions such as 
nc_inq_ncid().

\param ndimsp Pointer where number of dimensions will be
written. Ignored if NULL.

\returns ::NC_NOERR  No error.
\returns ::NC_EBADID Not a valid ID.

 */
int
nc_inq_ndims(int ncid, int *ndimsp)
{
    NC* ncp;
    int stat = NC_check_id(ncid, &ncp);
    if(stat != NC_NOERR) return stat;
    if(ndimsp == NULL) return NC_NOERR;
    return ncp->dispatch->inq(ncid,ndimsp,NULL,NULL,NULL);
}

/*!
Find the ID of the unlimited dimension.

This function finds the ID of the unlimited dimension. For
netCDF-4/HDF5 files (which may have more than one unlimited
dimension), the ID of the first unlimited dimesnion is returned. For
these files, nc_inq_unlimdims() will return all the unlimited dimension IDs.

\param ncid NetCDF or group ID, from a previous call to nc_open(),
nc_create(), nc_def_grp(), or associated inquiry functions such as 
nc_inq_ncid().

\param unlimdimidp Pointer where unlimited dimension ID will be
stored. If there is no unlimited dimension, -1 will be stored
here. Ignored if NULL.

\returns ::NC_NOERR  No error.
\returns ::NC_EBADID Not a valid ID.

 */
int
nc_inq_unlimdim(int ncid, int *unlimdimidp)
{
    NC* ncp;
    int stat = NC_check_id(ncid, &ncp);
    if(stat != NC_NOERR) return stat;
    return ncp->dispatch->inq_unlimdim(ncid,unlimdimidp);
}

/*!
Find out the name of a dimension.

\param ncid NetCDF or group ID, from a previous call to nc_open(),
nc_create(), nc_def_grp(), or associated inquiry functions such as 
nc_inq_ncid().

\param dimid Dimension ID, from a previous call to nc_inq_dimid() or
nc_def_dim().

\param name Returned dimension name. The caller must allocate space
for the returned name. The maximum possible length, in characters, of
a dimension name is given by the predefined constant
NC_MAX_NAME. (This doesn't include the null terminator, so declare
your array to be size NC_MAX_NAME+1). The returned character array
will be null-terminated. Ignored if NULL.

\returns ::NC_NOERR   No error.
\returns ::NC_EBADID  Not a valid ID.
\returns ::NC_EBADDIM Invalid dimension ID or name.

\section Example

Here is an example using nc_inq_dim() to determine the length of a
dimension named lat, and the name and current maximum length of the
unlimited dimension for an existing netCDF dataset named foo.nc:

\code
     #include <netcdf.h>
        ...
     int status, ncid, latid, recid;
     size_t latlength, recs;
     char recname[NC_MAX_NAME+1];
        ...
     status = nc_open("foo.nc", NC_NOWRITE, &ncid);
     if (status != NC_NOERR) handle_error(status);
     status = nc_inq_unlimdim(ncid, &recid); 
     if (status != NC_NOERR) handle_error(status);
        ...
     status = nc_inq_dimid(ncid, "lat", &latid);
     if (status != NC_NOERR) handle_error(status);
     status = nc_inq_dimlen(ncid, latid, &latlength);
     if (status != NC_NOERR) handle_error(status);

     status = nc_inq_dim(ncid, recid, recname, &recs);
     if (status != NC_NOERR) handle_error(status);
\endcode

 */
int
nc_inq_dimname(int ncid, int dimid, char *name)
{
    NC* ncp;
    int stat = NC_check_id(ncid, &ncp);
    if(stat != NC_NOERR) return stat;
    if(name == NULL) return NC_NOERR;
    return ncp->dispatch->inq_dim(ncid,dimid,name,NULL);
}

/*!
Find the length of a dimension.

The length for the unlimited dimension, if any, is the number of
records written so far.

\param ncid NetCDF or group ID, from a previous call to nc_open(),
nc_create(), nc_def_grp(), or associated inquiry functions such as 
nc_inq_ncid().

\param dimid Dimension ID, from a previous call to nc_inq_dimid() or
nc_def_dim().

\param lenp Pointer where the length will be stored.

\returns ::NC_NOERR   No error.
\returns ::NC_EBADID  Not a valid ID.
\returns ::NC_EBADDIM Invalid dimension ID or name.

\section Example

Here is an example using nc_inq_dim() to determine the length of a
dimension named lat, and the name and current maximum length of the
unlimited dimension for an existing netCDF dataset named foo.nc:

\code
     #include <netcdf.h>
        ...
     int status, ncid, latid, recid;
     size_t latlength, recs;
     char recname[NC_MAX_NAME+1];
        ...
     status = nc_open("foo.nc", NC_NOWRITE, &ncid);
     if (status != NC_NOERR) handle_error(status);
     status = nc_inq_unlimdim(ncid, &recid); 
     if (status != NC_NOERR) handle_error(status);
        ...
     status = nc_inq_dimid(ncid, "lat", &latid);  
     if (status != NC_NOERR) handle_error(status);
     status = nc_inq_dimlen(ncid, latid, &latlength); 
     if (status != NC_NOERR) handle_error(status);

     status = nc_inq_dim(ncid, recid, recname, &recs);
     if (status != NC_NOERR) handle_error(status);
\endcode
 */
int
nc_inq_dimlen(int ncid, int dimid, size_t *lenp)
{
    NC* ncp;
    int stat = NC_check_id(ncid, &ncp);
    if(stat != NC_NOERR) return stat;
    if(lenp == NULL) return NC_NOERR;
    return ncp->dispatch->inq_dim(ncid,dimid,NULL,lenp);
}
/**@}*/
