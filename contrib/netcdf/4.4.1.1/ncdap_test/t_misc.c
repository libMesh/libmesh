#include "config.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <netcdf.h>


#define URL1 "http://%s" /* test that no trailing / is ok */
static char url1[1024];

#ifdef DEBUG
static void
CHECK(int e, const char* msg)
{
    if(e == NC_NOERR) return;
    if(msg == NULL) msg = "Error";
    printf("%s: %s\n", msg, nc_strerror(e));
    exit(1);
}
#endif

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

    {
    char* evv = getenv("DTSTESTSERVER");
    if(evv == NULL)
	evv = "remotetest.unidata.ucar.edu";
    snprintf(url1,sizeof(url1),URL1,evv);
    }

    printf("Testing: Misc. Tests \n");
    retval = nc_open(url1, 0, &ncid);
    XFAIL(retval,"*** XFail : No trailing slash in url");
    retval = nc_close(ncid);
    return 0;
}
