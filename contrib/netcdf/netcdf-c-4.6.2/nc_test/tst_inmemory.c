/** \file \internal
Basic NC_INMEMORY API tests both for netcdf-3 and netcdf-4

Copyright 2011, UCAR/Unidata. See COPYRIGHT file for copying and
redistribution conditions.
*/

#undef DDBG

#include "config.h"
#include <stdio.h>
#include <stdlib.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#include "netcdf.h"
#include "netcdf_mem.h"
#include "ncbytes.h"
#include "nc_tests.h"
#include "err_macros.h"

#ifdef USE_NETCDF4
#include <hdf5.h>
extern int H5Eprint1(FILE * stream);
#endif

#define FLAGS4 (NC_INMEMORY|NC_NETCDF4|NC_CLOBBER)
#define FLAGS3 (NC_INMEMORY|NCCLOBBER)

#define NC_NETCDF3 0
#define MODIFIED 1

#define LARGE_SPACE (1<<18)

#define FILE3 "tst_inmemory3.nc"
#define FILE4 "tst_inmemory4.nc"
#define CREATE3 "tst_inmemory3_create.nc"
#define CREATE4 "tst_inmemory4_create.nc"
#define XFAIL "tst_xfail.nc"
#define MISC "tst_misc.nc"

#define CREATEFILE3 "tst_memcreate3.nc"
#define CREATEFILE4 "tst_memcreate4.nc"

/* Make no dimension larger than this */
#define MAXDIMLEN 100
#define NDIMS 2
#define UNLIM_LEN 2
#define DIM0_NAME "fun"
#define DIM1_NAME "money"
#define DIM1_LEN 8
#define ATT0_NAME "home"
#define ATT0_TEXT "earthship"
#define NVARS 3
#define VAR0_NAME "nightlife"
#define VAR1_NAME "time"
#define VAR2_NAME "taxi_distance"
#define VAR3_NAME "miles"

#define FLOATVAL ((float)42.22)

/*
CDL for created file:
netcdf tst_inmemory3 {
dimensions:
	fun = UNLIMITED ; // (2 currently)
	money = 8 ;
variables:
	int nightlife(fun, money) ;
	float time ;
	short taxi_distance(money) ;

// global attributes:
		:home = "earthship" ;
data:

 nightlife =
  0, 100, 200, 300, 400, 500, 600, 700,
  800, 900, 1000, 1100, 1200, 1300, 1400, 1500 ;

 time = 42.22 ;

 taxi_distance = 0, 1, 2, 3, 4, 5, 6, 7 ;
}
*/

#ifdef DDBG
#undef ERR
static void
fail(int line) {
    fflush(stdout);
    fprintf(stderr,"\nline=%d\n",line);
    fflush(stderr);
    exit(1);
}
#define ERR fail(__LINE__)
#endif

static int
check(int stat, const char* file, int line, int xfail)
{
/*    int pass = ((!xfail && stat == NC_NOERR) || (xfail && stat != NC_NOERR)); */
    fflush(stdout);
    if(!xfail) {
        if(stat != NC_NOERR) {
            fprintf(stderr,"***Fail: line: %d; status=%d %s\n",
	            line,stat,nc_strerror(stat));
        err++;
        } 
    } else { /*xfail*/
        if(stat == NC_NOERR) {
            fprintf(stderr,"***XFail Fail: line: %d; passed instead of failed\n",
	            line);
            err++;
        } else {
            fprintf(stderr,"\t***XFail: status=%d %s\n",
	            stat,nc_strerror(stat));
        } 
	stat = NC_NOERR; /* because xfail */
    }
    fflush(stderr);
    return stat;
}

#define CHECK(expr) {stat = check((expr),__FILE__,__LINE__,0); if(stat) return stat;}
#define XCHECK(expr) {stat = check((expr),__FILE__,__LINE__,1); if(stat) return stat;}

#define REPORT(xfail,expr) {if((xfail)) {XCHECK((expr));} else {CHECK((expr));}}

/**************************************************/

static void
removefile(const char* path)
{
    unlink(path);
}

static int
readfile(const char* path, NC_memio* memio)
{
    int status = NC_NOERR;
    FILE* f = NULL;
    size_t filesize = 0;
    size_t count = 0;
    char* memory = NULL;
    char* p = NULL;

    /* Open the file for reading */
#ifdef _MSC_VER
    f = fopen(path,"rb");
#else
    f = fopen(path,"r");
#endif
    if(f == NULL)
	{status = errno; goto done;}
    /* get current filesize */
    if(fseek(f,0,SEEK_END) < 0)
	{status = errno; goto done;}
    filesize = (size_t)ftell(f);
    /* allocate memory */
    memory = malloc((size_t)filesize);
    if(memory == NULL)
	{status = NC_ENOMEM; goto done;}
    /* move pointer back to beginning of file */
    rewind(f);
    count = filesize;
    p = memory;
    while(count > 0) {
        size_t actual;
        actual = fread(p,1,count,f);
	if(actual == 0 || ferror(f))
	    {status = NC_EIO; goto done;}	 
	count -= actual;
	p += actual;
    }
    if(memio) {
	memio->size = (size_t)filesize;
	memio->memory = memory;
    }    
done:
    if(status != NC_NOERR && memory != NULL)
	free(memory);
    if(f != NULL) fclose(f);
    return status;    
}


static int
writefile(const char* path, NC_memio* memio)
{
    int status = NC_NOERR;
    FILE* f = NULL;
    size_t count = 0;
    char* p = NULL;

    /* Open the file for writing */
#ifdef _MSC_VER
    f = fopen(path,"wb");
#else
    f = fopen(path,"w");
#endif
    if(f == NULL)
	{status = errno; goto done;}
    count = memio->size;
    p = memio->memory;
    while(count > 0) {
        size_t actual;
        actual = fwrite(p,1,count,f);
	if(actual == 0 || ferror(f))
	    {status = NC_EIO; goto done;}	 
	count -= actual;
	p += actual;
    }
done:
    if(f != NULL) fclose(f);
    return status;    
}


/* Duplicate an NC_memio instance; needed to avoid
   attempting to use memory that might have been realloc'd
   Allow the new memory to be larger than the src memory
*/
int
duplicatememory(NC_memio* src, NC_memio* target, size_t alloc)
{
    if(src == NULL || target == NULL || src->size == 0 || src->memory == NULL)
	return NC_EINVAL;
    *target = *src;
    if(alloc == 0) alloc = src->size;
    target->memory = malloc(alloc);
    if(target->memory == NULL)
	return NC_ENOMEM;
    memcpy(target->memory,src->memory,src->size);
    target->size = alloc;
    return NC_NOERR;
}

/*
Given an ncid of a created file, fill in the meta-data
and data as described by the above CDL. Do not close
the file.
*/

static int
define_metadata(int ncid)
{
    int stat = NC_NOERR;
    int dimid[NDIMS], varid0, varid1, varid2;
    short short_data[DIM1_LEN];
    size_t start[1] = {0};
    size_t count[1] = {DIM1_LEN};
    int dimprod = (UNLIM_LEN*DIM1_LEN);
    int i;
    float float_data;
    int nightdata[UNLIM_LEN*DIM1_LEN] ;

    /* Create data to write */
    float_data = FLOATVAL;

    for (i = 0; i < DIM1_LEN; i++)
        short_data[i] = i;

    for (i = 0; i < dimprod; i++)
        nightdata[i] = (100*i);

    CHECK(nc_put_att_text(ncid, NC_GLOBAL, ATT0_NAME,
			sizeof(ATT0_TEXT), ATT0_TEXT));

    CHECK(nc_def_dim(ncid, DIM0_NAME, NC_UNLIMITED, &dimid[0]));
    CHECK(nc_def_dim(ncid, DIM1_NAME, DIM1_LEN, &dimid[1]));

    CHECK(nc_def_var(ncid, VAR0_NAME, NC_INT, NDIMS, dimid, &varid0));
    CHECK(nc_def_var(ncid, VAR1_NAME, NC_FLOAT, 0, NULL, &varid1));
    CHECK(nc_def_var(ncid, VAR2_NAME, NC_SHORT, 1, &dimid[1], &varid2));

    CHECK(nc_enddef(ncid));

    CHECK(nc_put_vara_float(ncid, varid1, NULL, NULL, &float_data));
    CHECK(nc_put_vara_short(ncid, varid2, start, count, short_data));

    {
        size_t start[2] = {0,0};
        size_t count[2] = {2,DIM1_LEN};
        CHECK(nc_put_vara_int(ncid, varid0, start, count, nightdata));
    }

    return stat;
}


/*
Create our reference file as a real on-disk file
and read it back into memory
*/

static int
create_reference_file(const char* filename, int mode, NC_memio* filedata)
{
    int stat = NC_NOERR;
    int ncid;

    CHECK(nc_create(filename, mode|NC_CLOBBER, &ncid)); /* overwrite */
    CHECK(define_metadata(ncid));
    CHECK(nc_close(ncid));

    /* Read back the contents of the file into memory */
    if(filedata != NULL) {
	memset(filedata,0,sizeof(NC_memio));
	CHECK(readfile(filename,filedata));
    }
    return stat;
}

static int
modify_file(int ncid)
{
    int stat = NC_NOERR;
    size_t i;
    int varid3;
    int dimid[1];
    size_t len;
    int data[MAXDIMLEN];

    /* Get id of the unlimited dimension */
    if((stat=nc_inq_dimid(ncid, DIM0_NAME, dimid))) goto done;
    /* get current dim length */
    if((stat=nc_inq_dimlen(ncid, dimid[0], &len))) goto done;
    /* open file for new meta-data */
    if((stat=nc_redef(ncid))) goto done;
    /* Define a new variable */
    if((stat=nc_def_var(ncid, VAR3_NAME, NC_INT, 1, dimid, &varid3))) goto done;
    /* close metadata */
    if((stat=nc_enddef(ncid))) goto done;
    /* Write data to new variable */
    for(i=0;i<len;i++)
	data[i] = i;
    if((stat=nc_put_var_int(ncid,varid3,data))) goto done;
done:
    return stat;    
}

/* Verify the content of a file */
static int
verify_file(int ncid, int modified)
{
    int stat = NC_NOERR;
    int i;
    int dimid_in[NDIMS];
    int dimid[NDIMS];
    int ndims_in, nvars_in, natts_in, unlimdimid_in;
    char name_in[NC_MAX_NAME + 1], att0_in[NC_MAX_NAME + 1];
    nc_type type_in;
    size_t len_in;
    int varid[4];
    int nightdata_in[UNLIM_LEN*DIM1_LEN] ;
    float float_data_in;
    int milesdata_in[MAXDIMLEN];
    int dimprod = UNLIM_LEN * DIM1_LEN;
#ifdef USE_NETCDF4
    int tmp;
#endif
    
    CHECK(nc_inq(ncid, &ndims_in, &nvars_in, &natts_in, &unlimdimid_in));
    if (ndims_in != 2 || nvars_in != NVARS+modified || natts_in != 1 || unlimdimid_in != 0)
	CHECK(NC_EINVAL);

    /* Get all the dimids */
#ifdef USE_NETCDF4
    tmp = 0;
    CHECK((nc_inq_dimids(ncid,&tmp,dimid,1)));
    if(tmp != NDIMS) CHECK(NC_EINVAL);

    /* Get all the varids */
    tmp = 0;
    CHECK((nc_inq_varids(ncid,&tmp,varid)));
    if(tmp != (NVARS+modified)) CHECK(NC_EINVAL);
#else
    { /* Simulate nc_inq_varids and nc_inq_dimids */
	int j;
	int dimcnt = 0;
        int varcnt = 0;
	CHECK(nc_inq(ncid, &dimcnt, &varcnt, NULL, NULL));
	for(j=0;j<dimcnt;j++) dimid[j] = j;	
	for(j=0;j<varcnt;j++) varid[j] = j;	
    }
#endif

    CHECK(nc_get_att_text(ncid, NC_GLOBAL, ATT0_NAME, att0_in));
    att0_in[sizeof(ATT0_TEXT)] = '\0';
    if (strcmp(att0_in, ATT0_TEXT)) CHECK(NC_EINVAL);

    /* CHECK dimensions. */
    CHECK(nc_inq_dim(ncid, dimid[0], name_in, &len_in));
    if (strcmp(name_in, DIM0_NAME)) CHECK(NC_EINVAL);
    CHECK(nc_inq_dim(ncid, dimid[1], name_in, &len_in));
    if (strcmp(name_in, DIM1_NAME) || len_in != DIM1_LEN) CHECK(NC_EINVAL);

    /* CHECK variables. */
    CHECK(nc_inq_var(ncid, varid[0], name_in, &type_in, &ndims_in, dimid_in, &natts_in));
    if (strcmp(name_in, VAR0_NAME) || type_in != NC_INT || ndims_in != NDIMS ||
    dimid_in[0] != 0 || dimid_in[1] != 1 || natts_in != 0) CHECK(NC_EINVAL);
    CHECK(nc_inq_var(ncid, varid[1], name_in, &type_in, &ndims_in, dimid_in, &natts_in));
    if (strcmp(name_in, VAR1_NAME) || type_in != NC_FLOAT || ndims_in != 0 ||
    natts_in != 0) CHECK(NC_EINVAL);
    CHECK(nc_inq_var(ncid, varid[2], name_in, &type_in, &ndims_in, dimid_in, &natts_in));
    if (strcmp(name_in, VAR2_NAME) || type_in != NC_SHORT || ndims_in != 1 ||
        dimid_in[0] != 1 || natts_in != 0) CHECK(NC_EINVAL);

    CHECK(nc_get_var_int(ncid, varid[0], nightdata_in));
    for(i=0;i<dimprod;i++) {
	if(nightdata_in[i] != (100*i)) CHECK(NC_EINVAL);
    }

    CHECK(nc_get_vara_float(ncid, varid[1], NULL, NULL, &float_data_in));
    if (float_data_in != FLOATVAL) CHECK(NC_EINVAL);

    if(modified) {
	size_t unlimlen;
	CHECK(nc_inq_var(ncid, varid[3], name_in, &type_in, &ndims_in, dimid_in, &natts_in));
        if (strcmp(name_in, VAR3_NAME) || type_in != NC_INT || ndims_in != 1 ||
	    dimid_in[0] != 0 || natts_in != 0) CHECK(NC_EINVAL);
        CHECK(nc_inq_dimlen(ncid, dimid_in[0], &unlimlen));
        CHECK(nc_get_var_int(ncid, varid[3], milesdata_in));
	for(i=0;i<unlimlen;i++) {
	    if(milesdata_in[i] != i) CHECK(NC_EINVAL);
	}
    }

    return stat;
}

void
memiofree(NC_memio* memio)
{
    if(memio != NULL) {
	if(memio->memory != NULL)
	    free(memio->memory);
	memio->memory = NULL;
    }
}

static int
test_open(const char* path, NC_memio* filedata, int mode)
{
    int stat = NC_NOERR;
    NC_memio duplicate;
    NC_memio finaldata;
    int ncid;
    int xmode = mode; /* modified mode */

    finaldata.memory = NULL;
    finaldata.size = 0;
    finaldata.flags = 0;

    fprintf(stderr,"\n\t***Test open 1: nc_open_mem(): read-only\n");
    CHECK(duplicatememory(filedata,&duplicate,0));
    CHECK(nc_open_mem(path, xmode, duplicate.size, duplicate.memory, &ncid));
    CHECK(verify_file(ncid,!MODIFIED));
    CHECK(nc_close(ncid));
    memiofree(&duplicate);

    fprintf(stderr,"\n\t***Test open 2: nc_open_memio(): read-only\n");
    CHECK(duplicatememory(filedata,&duplicate,0));
    duplicate.flags = NC_MEMIO_LOCKED;
    CHECK(nc_open_memio(path, xmode, &duplicate, &ncid))
    CHECK(verify_file(ncid,!MODIFIED));
    CHECK(nc_close_memio(ncid,&finaldata));
    /* Published returned finaldata  */
    fprintf(stderr,"\tfinaldata: size=%lld memory=%p\n",(unsigned long long)finaldata.size,finaldata.memory);
    /* Verify that finaldata is same */
    if(finaldata.size != duplicate.size) CHECK(NC_EINVAL);
    if(finaldata.memory != duplicate.memory) CHECK(NC_EINVAL);
    memiofree(&finaldata);

    fprintf(stderr,"\n\t***Test open 3: nc_open_memio(): read-write, copy\n");
    xmode |= NC_WRITE; /* allow file to be modified */
    CHECK(duplicatememory(filedata,&duplicate,0));
    CHECK(nc_open_memio(path, xmode, &duplicate, &ncid))
    /* modify file */
    CHECK(modify_file(ncid));
    CHECK(verify_file(ncid,MODIFIED));
    CHECK(nc_close_memio(ncid,&finaldata));
    /* Published returned finaldata  */
    fprintf(stderr,"\tfinaldata: size=%lld memory=%p\n",(unsigned long long)finaldata.size,finaldata.memory);
    /* Verify that finaldata is same */
    if(finaldata.size < filedata->size) CHECK(NC_EINVAL);
    /* As a safeguard, the memory in duplicate should have been set to NULL*/
    memiofree(&finaldata);

    fprintf(stderr,"\n\t***Test open 4: nc_open_memio(): read-write, locked, extra space\n");
    /* Store the filedata in a memory chunk that leaves room for modification */
    CHECK(duplicatememory(filedata,&duplicate,LARGE_SPACE));
    /* Lock the duplicate memory */
    duplicate.flags |= NC_MEMIO_LOCKED;
    xmode |= NC_WRITE; /* allow file to be modified */
    CHECK(nc_open_memio(path, xmode, &duplicate, &ncid))
    /* modify file */
    CHECK(modify_file(ncid));
    CHECK(verify_file(ncid,MODIFIED));
    CHECK(nc_close_memio(ncid,&finaldata));
    /* Published returned finaldata  */
    fprintf(stderr,"\tfinaldata: size=%lld memory=%p\n",(unsigned long long)finaldata.size,finaldata.memory);
    /* Check returned finaldata:
       should have same memory but
       actual used final size should not exceed the original */
    if(finaldata.size > duplicate.size) CHECK(NC_EINVAL);
    if(finaldata.memory != duplicate.memory) CHECK(NC_EINVAL);
    memiofree(&finaldata);
    return stat;
}

static int
test_create(const char* path, int mode)
{
    int stat = NC_NOERR;
    NC_memio finaldata;
    int ncid;
    int xmode = mode;

    finaldata.memory = NULL;
    finaldata.size = 0;
    finaldata.flags = 0;

    fprintf(stderr,"\n\t***Test create 1: nc_create_memio(): no initialsize\n");
    CHECK(nc_create_mem(path, xmode, 0, &ncid))
    /* create file metadata */
    CHECK(define_metadata(ncid));
    CHECK(verify_file(ncid,!MODIFIED));
    CHECK(nc_close_memio(ncid,&finaldata));
    /* Published returned finaldata  */
    fprintf(stderr,"\tfinaldata: size=%lld memory=%p\n",(unsigned long long)finaldata.size,finaldata.memory);
    free(finaldata.memory);
    fprintf(stderr,"\n\t***Test create 2: nc_create_memio(): initialsize; save file\n");
    CHECK(nc_create_mem(path, xmode, LARGE_SPACE, &ncid))
    /* create file metadata */
    CHECK(define_metadata(ncid));
    CHECK(verify_file(ncid,!MODIFIED));
    CHECK(nc_close_memio(ncid,&finaldata));
    /* Published returned finaldata */
    fprintf(stderr,"\tfinaldata: size=%lld memory=%p\n",(unsigned long long)finaldata.size,finaldata.memory);
    /* Write out the final data as a .nc file */
    CHECK(writefile(path,&finaldata));
    if(finaldata.memory != NULL)
        free(finaldata.memory);
    finaldata.memory = NULL;
    return stat;
}

static int
test_misc(const char* path, int mode, NC_memio* filedata)
{
    int stat = NC_NOERR;
    int ncid;
    int xmode = mode;
    NC_memio duplicate;

    fprintf(stderr,"\n\t***Test misc 1: use nc_close on created inmemory file\n");
    CHECK(nc_create_mem(MISC, xmode, 0, &ncid))
    CHECK(nc_close(ncid));
    removefile(MISC);

    fprintf(stderr,"\n\t***Test misc 2: use nc_close on opened inmemory file\n");
    CHECK(duplicatememory(filedata,&duplicate,0));
    CHECK(nc_open_memio(path, xmode, &duplicate, &ncid))
    CHECK(verify_file(ncid,!MODIFIED));
    CHECK(nc_close(ncid));
    /* Do not free: nc_close will have done it memiofree(&duplicate); */
    removefile(MISC);

    return stat;
}

/* Test various edge conditions to ensure they fail correctly */
static int
test_xfail(const char* path, int mode, NC_memio* filedata)
{
    int stat = NC_NOERR;
    NC_memio duplicate = {0,NULL,0};
    int ncid;
    int xmode = mode; /* modified mode */

    fprintf(stderr,"\n\t***Test xfail 1: nc_open_mem(): write to read-only\n");
    CHECK(duplicatememory(filedata,&duplicate,0));
    CHECK(nc_open_mem(XFAIL, xmode, duplicate.size, duplicate.memory, &ncid));
    XCHECK(nc_redef(ncid));
    CHECK(nc_abort(ncid));
    memiofree(&duplicate);

    fprintf(stderr,"\n\t***Test xfail 2: nc_open_memio(): modify without overallocating\n");
    if((mode & NC_NETCDF4)) {
        fprintf(stderr,"\tSuppressed because of HDF5 library bug\n");
    } else {
      /* With HDF5 1.8.20, and possibly other versions,
         this tests causes a seg fault in the HDF5 Library.
         So until it is fixed, just leave well enough alone */
	NC_memio finaldata;
	memset(&finaldata,0,sizeof(finaldata));
	CHECK(duplicatememory(filedata,&duplicate,0));
	duplicate.flags = NC_MEMIO_LOCKED;
	xmode |= NC_WRITE;
	CHECK(nc_open_memio(XFAIL, xmode, &duplicate, &ncid))
	XCHECK(modify_file(ncid));    
	CHECK(nc_abort(ncid));
	memiofree(&finaldata);
	memiofree(&duplicate);
    }

    return stat;
}

int
main(int argc, char **argv)
{
    int stat = NC_NOERR;
    NC_memio filedata3;
#ifdef USE_NETCDF4
    NC_memio filedata4;
#endif

    fprintf(stderr,"\n*** Testing the inmemory API: netcdf-3.\n");
    CHECK(create_reference_file(FILE3,NC_NETCDF3,&filedata3)); /* netcdf-3 */
    CHECK(test_open(FILE3,&filedata3,NC_NETCDF3));
    CHECK(test_create(CREATE3,NC_NETCDF3));
    CHECK(test_misc(FILE3, NC_NETCDF3, &filedata3));
    CHECK(test_xfail(FILE3, NC_NETCDF3, &filedata3));
    memiofree(&filedata3);

#ifdef USE_NETCDF4
    fprintf(stderr,"\n*** Testing the inmemory API: netcdf-4.\n");
    CHECK(create_reference_file(FILE4,NC_NETCDF4,&filedata4));
    CHECK(test_open(FILE4,&filedata4,NC_NETCDF4));
    CHECK(test_create(CREATE4,NC_NETCDF4));
    CHECK(test_misc(FILE4,NC_NETCDF4, &filedata4));
    CHECK(test_xfail(FILE4, NC_NETCDF4, &filedata4));
    memiofree(&filedata4);
#endif

    SUMMARIZE_ERR;

    FINAL_RESULTS;

    return 0;
}
