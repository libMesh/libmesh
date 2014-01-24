/* This is part of the netCDF package. Copyright 2005 University
   Corporation for Atmospheric Research/Unidata See COPYRIGHT file for
   conditions of use. See www.unidata.ucar.edu for more info.

   Utility to rename a group 
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define LOGGING
#include "netcdf.h"

#undef RENAME_DEBUG

static void
usage(const char* msg)
{
    if(msg != NULL)
	fprintf(stderr,"%s\n",msg);
    fprintf(stderr,"usage: renamegroup <filename> <old group path name> <new name>\n");
    fflush(stderr);
    exit(0); /* do not cause error */
}

static void
check(int status)
{
    if(status == 0) return;
    fprintf(stderr,"%d: %s\n",status,nc_strerror(status));
    fflush(stderr);
    exit(1);
}

int
main(int argc, char **argv)
{
    int i,stat;
    int ncid, grpid;
    char* filename;
    char* oldname;
    char* newname;

    switch (argc) {
    case 0:
    case 1:
        usage("No arguments");
	break;
    case 2:
        usage("Too few arguments");
	break;
    case 3:
    default:
	filename = argv[1];
	oldname = argv[2];
	newname = argv[3];
	break;
    }

    if(strlen(filename) == 0)
	usage("bad filename argument");
    if(strlen(oldname) == 0)
	usage("bad old name argument");
    if(strlen(newname) == 0)
	usage("bad new name argument");

    stat = nc_open(filename,NC_WRITE,&ncid);
    check(stat);

#ifdef RENAME_DEBUG
    stat = nc_set_log_level(0);
    check(stat);
    nc_log_hdf5();
#endif

    stat = nc_inq_grp_full_ncid(ncid,oldname,&grpid);
    check(stat);

    stat = nc_rename_grp(grpid,newname);
    check(stat);

    stat = nc_close(ncid);
    check(stat);

    exit(0);
}

