#include "config.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <netcdf.h>


#define URL1 "http://remotetest.unidata.ucar.edu" /* test that no trailing / is ok */

static void
CHECK(int e, const char* msg)
{
    if(e == NC_NOERR) return;
    if(msg == NULL) msg = "Error";
    printf("%s: %s\n", msg, nc_strerror(e));
    exit(1);
}

static void
XFAIL(int e, const char* msg)
{
    if(e == NC_NOERR) return;
    if(msg == NULL) msg = "XFAIL";
    printf("%s: %s\n", msg, nc_strerror(e));
}

int
main()
{
    int ncid,retval;

    printf("Testing: Misc. Tests \n");
    retval = nc_open(URL1, 0, &ncid);
    XFAIL(retval,"*** XFail : No trailing slash in url");
    retval = nc_close(ncid);
    return 0;
}
