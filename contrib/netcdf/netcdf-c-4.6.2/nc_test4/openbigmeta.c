/*
Open a netcdf-4 file with horrendously large metadata.
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#include <netcdf.h>

#define FILE "bigmeta.nc"

int
main(int argc, char **argv)
{
    int ncid;
    time_t starttime, endtime;
    long long delta;

    starttime = 0;
    endtime = 0;
    time(&starttime);
    assert(nc_open(FILE,NC_NETCDF4,&ncid) == NC_NOERR);
    time(&endtime);
    assert(nc_close(ncid) == NC_NOERR);

    /* Compute the delta 1 second resolution is fine for this */
    delta = (long long)(endtime - starttime);
    printf("open delta=%lld\n",delta);
    return 0;
}
