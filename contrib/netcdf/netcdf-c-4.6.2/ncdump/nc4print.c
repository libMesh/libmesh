/*********************************************************************
 *   Copyright 2016, UCAR/Unidata
 *   See netcdf/COPYRIGHT file for copying and redistribution conditions.
 *********************************************************************/

/**
This provides a simple netcdf-4 metadata -> xml printer.
Primarily for use in debugging, but could be adapted to
create other tools.
*/

/**************************************************/

#include "config.h"
#include <stdlib.h>
#include <stdio.h>
#include "netcdf.h"
#include "ncbytes.h"

extern int NC4print(NCbytes* buf, int ncid);

int
main(int argc, char** argv)
{
    int i;
    int ret = NC_NOERR;

    if(argc == 1) {
	fprintf(stderr,"usage: nc4printer <file> <file>...\n");
	exit(1);
    }
    for(i=1;i<argc;i++) {
	int ncid;
	char* filename;
        NCbytes* buf;

	filename = argv[i];
	buf = ncbytesnew();

	if((ret = nc_open(filename,NC_NETCDF4,&ncid))) goto fail;
		
	ret = NC4print(buf,ncid);
	ncbytesnull(buf);
	fprintf(stderr,"========== %s ==========\n",filename);
	fprintf(stderr,"%s\n",ncbytescontents(buf));
	ncbytesfree(buf);
    }
    exit(0);

fail:
    fprintf(stderr,"***Fail: (%d) %s\n",ret,nc_strerror(ret));
    exit(1);
}
