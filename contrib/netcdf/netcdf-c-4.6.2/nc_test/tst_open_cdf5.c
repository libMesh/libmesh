#include "config.h"
#include <stdio.h>
#include <stdlib.h>
#include <netcdf.h>

#define FILE_NAME "bad_cdf5_begin.nc"

int main(int argc, char *argv[])
{
    char *fname=FILE_NAME;
    int err, nerrs=0, ncid;

    if (argc == 2) fname = argv[1];

    err = nc_open(fname, NC_NOWRITE, &ncid);
    if (err != NC_ENOTNC) {
        printf("Error: nc_open() expect NC_ENOTNC but got (%s)\n",
               nc_strerror(err));
        nerrs++;
    }
    else if (err == NC_NOERR) /* close file */
        nc_close(ncid);

    return (nerrs > 0);
}
