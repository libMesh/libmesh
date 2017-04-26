/*
 * Test contributed in support of netCDF issue
 * https://github.com/Unidata/netcdf-c/issues/160
 */

#include "netcdf.h"

#include <stdio.h>

int main(int argc, char* argv[])
{
    int ret;
    int ncid;
    int dim1id, dim2id;
    int var1id, var2id;
    size_t start = 0;
    size_t count = 5;
    double vals[] = { 1.0, 2.0, 3.0, 4.0, 5.0 };

    if ((ret = nc_create("tst_put_vars_two_unlim_dim.nc", NC_NETCDF4 | NC_CLOBBER, &ncid))) {
        printf("nc_create(...): error code = %d\n", ret);
        return -1;
    }

    if ((ret = nc_def_dim(ncid, "dim1", NC_UNLIMITED, &dim1id))) {
        printf("nc_def_dim(...\"dim1\"...): error code = %d\n", ret);
        nc_close(ncid);
        return -1;
    }

    if ((ret = nc_def_dim(ncid, "dim2", NC_UNLIMITED, &dim2id))) {
        printf("nc_def_dim(...\"dim1\"...): error code = %d\n", ret);
        nc_close(ncid);
        return -1;
    }

    if ((ret = nc_def_var(ncid, "var1", NC_DOUBLE, 1, &dim1id, &var1id))) {
        printf("nc_def_var(...\"var1\"...): error code = %d\n", ret);
        nc_close(ncid);
        return -1;
    }

    if ((ret = nc_def_var(ncid, "var2", NC_DOUBLE, 1, &dim2id, &var2id))) {
        printf("nc_def_var(...\"var2\"...): error code = %d\n", ret);
        nc_close(ncid);
        return -1;
    }

    if ((ret = nc_put_vars_double(ncid, var1id, &start, &count, NULL, &vals[0])))  {
        printf("nc_put_var_double(...var1id...): error code = %d\n", ret);
        nc_close(ncid);
        return -1;
    }

    if ((ret = nc_put_vars_double(ncid, var2id, &start, &count, NULL, &vals[0])))  {
        printf("nc_put_var_double(...var2id...): error code = %d\n", ret);
        nc_close(ncid);
        return -1;
    }

    if ((ret = nc_close(ncid))) {
        printf("nc_close(...): error code = %d\n", ret);
        return -1;
    }

    return 0;
}
