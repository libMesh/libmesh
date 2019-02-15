/*********************************************************************
 *   Copyright 2016, UCAR/Unidata
 *   See netcdf/COPYRIGHT file for copying and redistribution conditions.
 *********************************************************************/

/**
Test the netcdf-4 data building process.
*/

#include <stdlib.h>
#include <stdio.h>
#include "netcdf.h"

static void
fail(int code)
{
    if(code != NC_NOERR)
	fprintf(stderr,"***Fail: %s\n",nc_strerror(code));
    exit((code==NC_NOERR?EXIT_SUCCESS:EXIT_FAILURE));
}

int
main(int argc, char** argv)
{
    int ret = NC_NOERR;
    char url[4096];
    int ncid;

    /* Skip cmd name */
    argc++;
    argv++;

    if(argc < 2) {
	fprintf(stderr, "too few arguments: t_dmrdata.c <infile> <outfile>\n");
	fail(NC_NOERR);
    }

    /* build the url */
    snprintf(url,sizeof(url),"file://%s#dap4&debug=copy&substratename=%s",argv[0],argv[1]);

#ifdef DEBUG
    fprintf(stderr,"t_dmrbuild %s -> %s\n",url,outfile);
#endif
  
    /* Use the open/close mechanism */
    if((ret=nc_open(url,NC_NETCDF4,&ncid))) goto done;
    if((ret=nc_close(ncid))) goto done;

done:
#ifdef DEBUG
    fprintf(stderr,"code=%d %s\n",ret,nc_strerror(ret));
#endif
    return (ret ? 1 : 0);
}
