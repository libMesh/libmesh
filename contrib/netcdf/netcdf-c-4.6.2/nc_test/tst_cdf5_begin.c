#include <stdio.h>
#include <stdlib.h>
#include <netcdf.h>

/* When using NetCDF 4.4.1 ad prior to create a CDF-5 file and defining a small
 * variable after a big variable (> 2^32-3 bytes), the file starting offset of
 * the small variable (and all variables defined after the big variable) is
 * calculated incorrectly. This test program detects this bug by checking the
 * contents of the possible overlaps between the two variables.
 */

#define ERR {if(err!=NC_NOERR){printf("Error at line %d in %s: %s\n", __LINE__,__FILE__, nc_strerror(err));nerrs++;}}

#define FILE_NAME "tst_cdf5_begin.nc"

int main(int argc, char *argv[])
{
    int i, err, nerrs=0, ncid, dimid[2], varid[2];
    short buf[10];
    size_t start, count;
  
    err = nc_create(FILE_NAME, NC_CLOBBER|NC_64BIT_DATA, &ncid); ERR;
    err = nc_def_dim(ncid, "dim0", NC_MAX_UINT, &dimid[0]); ERR
    err = nc_def_dim(ncid, "dim1", 10,          &dimid[1]); ERR

    /* define one small variable after one big variable */
    err = nc_def_var(ncid, "var_big",   NC_SHORT, 1, &dimid[0], &varid[0]); ERR
    err = nc_def_var(ncid, "var_small", NC_SHORT, 1, &dimid[1], &varid[1]); ERR
    err = nc_set_fill(ncid, NC_NOFILL, NULL); ERR
    err = nc_enddef(ncid); ERR

    /* write to var_big in location overlapping with var_small when using
     * netCDF 4.4.x or prior */
    start = NC_MAX_UINT/sizeof(short);
    count = 10;
    for (i=0; i<10; i++) buf[i] = i;
    err = nc_put_vara_short(ncid, varid[0], &start, &count, buf); ERR

    /* write var_small */
    for (i=0; i<10; i++) buf[i] = -1;
    err = nc_put_var_short(ncid, varid[1], buf); ERR

    /* read back var_big and check contents */
    for (i=0; i<10; i++) buf[i] = -1;
    err = nc_get_vara_short(ncid, varid[0], &start, &count,buf); ERR
    for (i=0; i<10; i++) {
        if (buf[i] != i) {
            printf("Error at buf[%d] expect %d but got %hd\n",i,i,buf[i]);
            nerrs++;
        }
    }
    err = nc_close(ncid); ERR

    return (nerrs > 0);
}

