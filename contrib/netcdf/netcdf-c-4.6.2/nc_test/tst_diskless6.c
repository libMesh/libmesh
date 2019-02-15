#include "config.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#ifdef HAVE_TIME_H
#include <time.h>
#endif
#ifdef HAVE_SYS_STAT_H
#include <sys/stat.h>
#endif
#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>
#endif
#ifdef HAVE_UNISTD_H
#include "unistd.h"
#endif
#ifdef _WIN32
#include <Windows.h>
#include <io.h>
#endif

#include "netcdf.h"

#define CLEANUP

#ifndef NC_NETCDF3
#define NC_NETCDF3 0
#endif

#define FILE3D "file3d.nc"
#define FILE3DP "file3dp.nc"
#define FILE4D "file4d.nc"
#define FILE4DP "file4dp.nc"

/* Mnemonics */
#define RDONLY 1
#define RDWRITE 0

#ifdef _WIN32
#define RDONLYMODE (_S_IREAD)
#define RDWRMODE (_S_IREAD|_S_IWRITE)
#define CHMOD _chmod
#define SLEEP(x) Sleep((x)*1000)
#else
#define RDONLYMODE (S_IRUSR|S_IRGRP|S_IROTH)
#define RDWRMODE (RDONLY | S_IWUSR)
#define CHMOD chmod
#define SLEEP(x) sleep(x)
#endif

typedef enum OC { OPEN, CLOSE} OC;

static int lineno = 0;

static char*
smode(int mode)
{
    static char ms[8192];
    ms[0] = '\0';
    if(mode & NC_NETCDF4)
	strcat(ms,"NC_NETCDF4");
    else
	strcat(ms,"NC_NETCDF3");
    if(mode & NC_DISKLESS)
	strcat(ms,"|NC_DISKLESS");
    if(mode & NC_WRITE)
	strcat(ms,"|NC_WRITE");
    if(mode & NC_NOCLOBBER)
	strcat(ms,"|NC_NOCLOBBER");
    if(mode & NC_INMEMORY)
	strcat(ms,"|NC_INMEMORY");
    if(mode & NC_PERSIST)
	strcat(ms,"|NC_PERSIST");
    if(mode & NC_MMAP)
	strcat(ms,"|NC_MMAP");
    return ms;
}

/* Return 1 if file was changed else 0 */
static time_t
getmodified(const char* path)
{
    struct stat attr;
    stat(path, &attr);
    return attr.st_mtime;
}

static void
changeaccess(int rdonly)
{
   const char* p3 = FILE3DP; /* Keep Visual Studio happy */
   const char* p4 = FILE4DP;
   int mode = (rdonly?RDONLYMODE:RDWRMODE);
   (void)CHMOD(p3,mode);
   (void)CHMOD(p4,mode);
}

static void
cleanup()
{
    changeaccess(RDWRITE);
    /* cleanup */
    (void)unlink(FILE3DP);
    (void)unlink(FILE4DP);
}

static void
fail(int ret)
{
    if(ret != NC_NOERR) {
        fprintf(stderr,"*** Fail: line: %d: (%d) %s\n", lineno, ret, nc_strerror(ret));
        fflush(stderr);
#ifdef CLEANUP
	cleanup();
#endif
        exit(1);
    }
}

void
exists(const char* file)
{
    FILE* f = fopen(file,"r");
    if(f == NULL) fail(NC_EPERM);
    fclose(f);
}

void
notexists(const char* file)
{
    FILE* f = fopen(file,"r");
    if(f != NULL) {fclose(f); fail(NC_EEXIST);}
}

#define TESTCREATE(file,mode,exist) {lineno=__LINE__; testcreate(file,mode,exist);}
#define TESTOPEN(file,mode,rw) {lineno=__LINE__; testopen(file,mode,rw);}

static void
testcreate(const char* file, int mode, int mustexist)
{
    int ret = NC_NOERR;
    int ncid, dimid;

    printf("test: file=%s mode=%s\n",file,smode(mode)); fflush(stdout);

    if((ret = nc_create(file,mode,&ncid))) fail(ret);
    if((ret = nc_def_dim(ncid,"dim",5,&dimid))) fail(ret);
    if((ret = nc_close(ncid))) fail(ret);
    if(mustexist)
        exists(file);
    else
        notexists(file);
}

static void
testopen(const char* file, int mode, int rdwrite)
{
    int ret = NC_NOERR;
    int ncid, dimid;
    size_t len;
    time_t time1, time2;

    printf("test: file=%s mode=%s\n",file,smode(mode)); fflush(stdout);

    time1 = getmodified(file);
    SLEEP(1); /* Ensure that if modified, it will be reported */
    if((ret = nc_open(file,mode,&ncid))) fail(ret);
    if(rdwrite) {
        /* Modify */
        if((ret=nc_redef(ncid))) fail(ret);
	/* See if dim2 is already defined */
        ret = nc_inq_dimid(ncid,"dim2",&dimid);
	if(ret == NC_NOERR) {/* defined */
            if((ret = nc_def_dim(ncid,"dim3",3,&dimid))) fail(ret);
	} else if(ret == NC_EBADDIM) { /* not defined */
            if((ret = nc_def_dim(ncid,"dim2",2,&dimid))) fail(ret);
	} else
	   fail(ret);
        if((ret=nc_enddef(ncid))) fail(ret);
    }
    if((ret = nc_inq_dimid(ncid,"dim",&dimid))) fail(ret);
    if((ret = nc_inq_dimlen(ncid,dimid,&len))) fail(ret);
    ret = nc_close(ncid);
    if(ret != NC_NOERR && rdwrite) fail(ret);
    else  if(ret == NC_NOERR && !rdwrite) fail(NC_EINVAL);
    time2 = getmodified(file);
    if(!rdwrite) {
	if(time2 != time1) {
	    fprintf(stderr,"file modified: %s\n",file);
	    fail(NC_EACCESS);
	}
    }
}

int
main()
{
    changeaccess(RDWRITE);
    cleanup();

    /* Test various mode combinations */

    /* Create and persist some files using diskless */
    printf("*** Test create\n"); fflush(stdout);
    TESTCREATE(FILE3D,NC_NETCDF3|NC_DISKLESS,0); /* Not persisted */
    TESTCREATE(FILE3DP,NC_NETCDF3|NC_DISKLESS|NC_PERSIST,1); /* persisted */
    TESTCREATE(FILE4D,NC_NETCDF4|NC_DISKLESS,0); /* Not persisted */
    TESTCREATE(FILE4DP,NC_NETCDF4|NC_DISKLESS|NC_PERSIST,1); /* persisted */

    /* Of the above persisted files, re-read and modify and re-persist */
    printf("*** Test open + modify + rdwrite\n"); fflush(stdout);
    TESTOPEN(FILE3DP,NC_NETCDF3|NC_DISKLESS|NC_PERSIST|NC_WRITE,1); /* re-persist */
    TESTOPEN(FILE4DP,NC_NETCDF4|NC_DISKLESS|NC_PERSIST|NC_WRITE,1); /* re-persist */

    /* Of the above modified files, re-read and modify but do not re-persist */
    /* Test open + modify + rdonly; requires NC_DISKLESS*/
#if 0
    /* Fails with hdf5, file must be writeable even if not persisted */ 
    changeaccess(RDONLY); /* prevent re-persist */
#endif
    printf("*** Testopen modify + rdonly\n"); fflush(stdout);
    TESTOPEN(FILE3DP,NC_NETCDF3|NC_DISKLESS|NC_WRITE,1);
    TESTOPEN(FILE4DP,NC_NETCDF4|NC_DISKLESS|NC_WRITE,1);

#ifdef CLEANUP
    cleanup();
#endif
    return 0;
}
