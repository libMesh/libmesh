#include <stdio.h>
#include <netcdf.h>

#define ERR {if(err!=NC_NOERR){printf("Error at line %d in %s: %s\n", __LINE__,__FILE__, nc_strerror(err));nerrs++;}}

#define EXP_ERR(exp,err) { \
    if (err != exp) { \
        nerrs++; \
        printf("Error at line %d in %s: expecting %s but got %d\n", \
        __LINE__,__FILE__,#exp, err); \
    } \
}

int main(int argc, char *argv[])
{
    int i, err, nerrs=0, ncid, dimid[NC_MAX_VAR_DIMS+2], varid;

    err = nc_create("tst_max_var_dims.nc", NC_CLOBBER, &ncid); ERR;
    err = nc_def_dim(ncid, "dim0", NC_UNLIMITED, &dimid[0]); ERR;
    err = nc_def_dim(ncid, "dim1", 1, &dimid[1]); ERR;

    for (i=2; i<NC_MAX_VAR_DIMS+2; i++) dimid[i] = dimid[1];

    err = nc_def_var(ncid, "v0", NC_INT, NC_MAX_VAR_DIMS+1, &dimid[0], &varid);
    EXP_ERR(NC_EMAXDIMS,err)

    err = nc_def_var(ncid, "v1", NC_INT, NC_MAX_VAR_DIMS+1, &dimid[1], &varid);
    EXP_ERR(NC_EMAXDIMS,err)

    err = nc_set_fill(ncid, NC_NOFILL, NULL); ERR
    err = nc_close(ncid); ERR
    return (nerrs > 0);
}

