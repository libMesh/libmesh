#include <stdio.h>
#include <netcdf.h>

#define CHECK_ERR { \
    if (err != NC_NOERR) { \
        nerrs++; \
        printf("Error at line %d in %s: (%s)\n", \
        __LINE__,__FILE__,nc_strerror(err)); \
    } \
}

#define EXP_ERR(exp) { \
    if (err != exp) { \
        nerrs++; \
        printf("Error at line %d in %s: expecting %s but got %s\n", \
        __LINE__,__FILE__,nc_strerror(exp), nc_strerror(err)); \
    } \
}

int main(int argc, char** argv)
{
    char *filename="tst_err_enddef.nc";
    int err, nerrs=0, ncid, cmode, varid, dimid[3];

    if (argc == 2) filename = argv[1];
    printf("*** TESTING error code returned from nc__enddef and nc_close ");

    cmode = NC_CLOBBER;
    err = nc_create(filename, cmode, &ncid); CHECK_ERR
    err = nc_set_fill(ncid, NC_NOFILL, NULL); CHECK_ERR

    err = nc_def_dim(ncid, "X", 5,      &dimid[0]); CHECK_ERR
    err = nc_def_dim(ncid, "YY", 32000, &dimid[1]); CHECK_ERR
    err = nc_def_dim(ncid, "XX", 32000, &dimid[2]); CHECK_ERR

    err = nc_def_var(ncid, "var", NC_INT, 1, dimid, &varid); CHECK_ERR
    err = nc_def_var(ncid, "var_big", NC_FLOAT, 2, dimid+1, &varid); CHECK_ERR

    /* make the file header size larger than 2 GiB */
    err = nc__enddef(ncid, 2147483648LL, 1, 1, 1);
    EXP_ERR(NC_EVARSIZE)

    /* the above error keeps the program in define mode, thus close will
     * call enddef again, but this time no error is expected
     */
    err = nc_close(ncid); CHECK_ERR

    if (nerrs) printf(".... failed with %d errors\n",nerrs);
    else       printf(".... pass\n");

    return (nerrs > 0);
}
