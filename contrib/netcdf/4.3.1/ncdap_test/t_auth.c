#include "config.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <netcdf.h>

#define BASICAUTHURL "http://tiggeUser:tigge@thredds-test.ucar.edu/thredds/dodsC/restrict/testData.nc"

static void
CHECK(int e, const char* msg)
{
    if(e == NC_NOERR) return;
    if(msg == NULL) msg = "Error";
    printf("%s: %s\n", msg, nc_strerror(e));
    exit(1);
}

int
main()
{
    int ncid,retval;

    printf("Testing: Http Basic Authorization\n");
    retval = nc_open(BASICAUTHURL, 0, &ncid);
    CHECK(retval,"*** Fail: Http Basic Authorization");
    retval = nc_close(ncid);
    printf("*** PASS: Http Basic Authorization\n");
    return 0;
}
