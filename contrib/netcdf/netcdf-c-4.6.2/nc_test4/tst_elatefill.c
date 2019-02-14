/* This is part of the netCDF package. Copyright 2017 University
   Corporation for Atmospheric Research/Unidata See COPYRIGHT file for
   conditions of use. See www.unidata.ucar.edu for more info.

   Test proper elatefill return when fillvalue is assigned outside of
   the initial define.

   Contributed by wkliao, see the following for more information:

   * https://github.com/Unidata/netcdf-c/issues/384
   * https://github.com/Unidata/netcdf-c/pull/387
   * https://github.com/Unidata/netcdf-c/issues/390
*/

#include "config.h"
#include <nc_tests.h>
#include <stdio.h>
#include <netcdf.h>

#define FILE_NAME "tst_elatefill.nc"

#define ERR_CHK {if(err!=NC_NOERR)printf("Error at line %d: %s\n",__LINE__,nc_strerror(err));}

int
main(int argc, char **argv)
{
    int ncid, dimid, varid, err;
    int fillv;

    err = nc_create(FILE_NAME, NC_NETCDF4, &ncid); ERR_CHK;
    err = nc_def_dim(ncid, "dim", 10, &dimid); ERR_CHK;
    err = nc_def_var(ncid, "var", NC_INT, 1, &dimid, &varid); ERR_CHK;
    err = nc_enddef(ncid); ERR_CHK;

    err = nc_redef(ncid); ERR_CHK;

    /* try put attribute _FillValue and expect NC_ELATEFILL */
    fillv = 9;
    err = nc_put_att_int(ncid, varid, _FillValue, NC_INT, 1, &fillv);
    if (err != NC_ELATEFILL)
        printf("line %d expecting NC_ELATEFILL but got %d\n",__LINE__,err);
    err = nc_close(ncid); ERR_CHK;
    return 0;
}
