#include "config.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "netcdf.h"

/* Embedded user:pwd */
static char* URL1 = 
"http://tiggeUser:tigge@remotetest.unidata.ucar.edu/thredds/dodsC/restrict/testData.nc";
/* user:pwd from .dodsrc*/
static char* URL2 = 
"http://remotetest.unidata.ucar.edu/thredds/dodsC/restrict/testData.nc";

/* .dodsrc file */
static char* CONTENT = "HTTP.CREDENTIALS.USER=tiggeUser\nHTTP.CREDENTIALS.PASSWORD=tigge\n";

static void
CHECK(int e, const char* msg)
{
    if(e == NC_NOERR) return;
    if(msg == NULL) msg = "Error";
    printf("%s: %s\n", msg, nc_strerror(e));
    exit(1);
}


int
main(int argc, char** argv)
{
    int ncid,retval,pass;
    char** url;
    FILE* dodsrc;
    pass = 1; /* assume success */

    dodsrc = fopen(".dodsrc","w");
    if(dodsrc == NULL) {
        fprintf(stderr,"Cannot create .dodsrc\n");
        exit(1);
    }    
    fprintf(dodsrc,CONTENT);
    fclose(dodsrc);

    printf("Testing: Http Basic Authorization\n\n");
    if(1) {
        printf("Embedded user:pwd: %s\n",URL1);
        retval = nc_open(URL1, 0, &ncid);
        if(retval != NC_NOERR) {
            pass = 0;
            printf("*** FAIL: Embedded user:pwd\n");
        } else {
            printf("*** PASS: Embedded user:pwd\n");
	    retval = nc_close(ncid);
	}
        fflush(stdout);
    }

    if(1) {
        printf(".dodsrc user:pwd: %s\n",URL1);

        retval = nc_open(URL2, 0, &ncid);
        if(retval != NC_NOERR) {
            pass = 0;
            printf("*** FAIL: .dodsrc user:pwd\n");
        } else {
	    retval = nc_close(ncid);
            printf("*** PASS: .dodsrc user:pwd\n");
        }
        fflush(stdout);
    }
    unlink(".dodsrc"); /* delete the file */
    return !pass;
}
