/*
  Copyright 2008, UCAR/Unidata
  See COPYRIGHT file for copying and redistribution conditions.
*/

#include "config.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <hdf5.h>
#include "netcdf.h"

#undef DEBUG

/* The C standard apparently defines all floating point constants as double;
   we rely on that in this code.
*/
#define DBLVAL 12345678.12345678

#define TEST_ID 32768

#define MAXERRS 8

#define MAXPARAMS 32

#define NPARAMS 14

static unsigned int baseline[NPARAMS];

#define MAXDIMS 8

#define DEFAULTACTUALDIMS 4
#define DEFAULTDIMSIZE 4
#define DEFAULTCHUNKSIZE 4

#define TESTFILE "testmisc.nc"

#define spec "32768, -17b, 23ub, -25S, 27US, 77, 93U, 789f, 12345678.12345678d, -9223372036854775807L, 18446744073709551615UL"

static size_t dimsize = DEFAULTDIMSIZE;
static size_t chunksize = DEFAULTCHUNKSIZE;
static size_t actualdims = DEFAULTACTUALDIMS;
static size_t pattern[MAXDIMS];

static size_t totalproduct = 1; /* x-product over max dims */
static size_t actualproduct = 1; /* x-product over actualdims */
static size_t chunkproduct = 1; /* x-product over actual chunks */

static size_t dims[MAXDIMS];
static size_t chunks[MAXDIMS];

static int nerrs = 0;

static int ncid, varid;
static int dimids[MAXDIMS];
static size_t odom[MAXDIMS];
static float* array = NULL;
static float* expected = NULL;

static unsigned int filterid = 0;
static size_t nparams = 0;
static unsigned int params[MAXPARAMS];

/* Forward */
static int test_test1(void);
static void init(int argc, char** argv);
static void reset(void);
static void odom_reset(void);
static int odom_more(void);
static int odom_next(void);
static int odom_offset(void);
static float expectedvalue(void);
static void verifyparams(void);

#define ERRR do { \
fflush(stdout); /* Make sure our stdout is synced with stderr. */ \
fprintf(stderr, "Sorry! Unexpected result, %s, line: %d\n", \
        __FILE__, __LINE__);                                \
nerrs++;\
} while (0)

static int
check(int err,int line)
{
    if(err != NC_NOERR) {
        fprintf(stderr,"fail (%d): %s\n",line,nc_strerror(err));
    }
    return NC_NOERR;
}

static void
report(const char* msg, int lineno)
{
    fprintf(stderr,"fail: line=%d %s\n",lineno,msg);
    exit(1);
}

#define CHECK(x) check(x,__LINE__)
#define REPORT(x) report(x,__LINE__)

static int
verifychunks(void)
{
    int i;
    int store = -1;
    size_t chunksizes[MAXDIMS];
    memset(chunksizes,0,sizeof(chunksizes));
    CHECK(nc_inq_var_chunking(ncid, varid, &store, chunksizes));
    if(store != NC_CHUNKED) {
        fprintf(stderr,"bad chunk store\n");
        return 0;
    }
    for(i=0;i<actualdims;i++) {
        if(chunksizes[i] != chunks[i]) {
            fprintf(stderr,"bad chunk size: %d\n",i);
            return 0;
        }
    }
    return 1;
}

static int
create(void)
{
    int i;

    /* Create a file with one big variable. */
    CHECK(nc_create(TESTFILE, NC_NETCDF4|NC_CLOBBER, &ncid));
    CHECK(nc_set_fill(ncid, NC_NOFILL, NULL));
    for(i=0;i<actualdims;i++) {
        char dimname[1024];
        snprintf(dimname,sizeof(dimname),"dim%d",i);
        CHECK(nc_def_dim(ncid, dimname, dims[i], &dimids[i]));
    }
    CHECK(nc_def_var(ncid, "var", NC_FLOAT, actualdims, dimids, &varid));
    return NC_NOERR;
}

static void
setvarfilter(void)
{
    CHECK(nc_def_var_filter(ncid,varid,TEST_ID,NPARAMS,baseline));
    verifyparams();
}

static void
verifyparams(void)
{
    int i;
    CHECK(nc_inq_var_filter(ncid,varid,&filterid,&nparams,params));
    if(filterid != TEST_ID) REPORT("id mismatch");
    if(nparams != NPARAMS) REPORT("nparams mismatch");
    for(i=0;i<nparams;i++) {
        if(params[i] != baseline[i])
            REPORT("param mismatch");
    }
}

static int
openfile(void)
{
    unsigned int* params;

    /* Open the file and check it. */
    CHECK(nc_open(TESTFILE, NC_NOWRITE, &ncid));
    CHECK(nc_inq_varid(ncid, "var", &varid));

    /* Check the compression algorithm */
    CHECK(nc_inq_var_filter(ncid,varid,&filterid,&nparams,NULL));
    if(nparams > 0) {
        params = (unsigned int*)malloc(sizeof(unsigned int)*nparams);
        if(params == NULL)
            return NC_ENOMEM;
        CHECK(nc_inq_var_filter(ncid,varid,&filterid,&nparams,params));
    }
    if(filterid != TEST_ID) {
        fprintf(stderr,"open: test id mismatch: %d\n",filterid);
        return NC_EFILTER;
    }
    if(nparams != NPARAMS) {
	size_t i;
	unsigned int inqparams[MAXPARAMS];
        fprintf(stderr,"nparams  mismatch\n");
        for(nerrs=0,i=0;i<nparams;i++) {
            if(inqparams[i] != baseline[i]) {
                fprintf(stderr,"open: testparam mismatch: %ld\n",(unsigned long)i);
		nerrs++;
	    }
	}
    }
    if(nerrs > 0) return NC_EFILTER; 

    /* Verify chunking */
    if(!verifychunks())
        return 0;
    fflush(stderr);
    return 1;
}

static int
setchunking(void)
{
    int store;

    store = NC_CHUNKED;
    CHECK(nc_def_var_chunking(ncid,varid,store,chunks));
    if(!verifychunks())
        return NC_EINVAL;
    return NC_NOERR;
}

static void
fill(void)
{
   odom_reset();
   if(1) {
        int i;
        if(actualproduct <= 1) abort();
        for(i=0;i<actualproduct;i++)
            expected[i] = (float)i;
   } else {
       while(odom_more()) {
            int offset = odom_offset();
            float expect = expectedvalue();
            expected[offset] = expect;
            odom_next();
        }
   }
}


static int
compare(void)
{
    int errs = 0;
    fprintf(stderr,"data comparison: |array|=%ld\n",(unsigned long)actualproduct);
    if(1)
    {
        int i;
        for(i=0;i<actualproduct;i++) {
            if(expected[i] != array[i]) {
                fprintf(stderr,"data mismatch: array[%d]=%f expected[%d]=%f\n",
                            i,array[i],i,expected[i]);
                errs++;
                if(errs >= MAXERRS)
                    break;
            }
        }
   } else
   {
       odom_reset();
       while(odom_more()) {
            int offset = odom_offset();
            float expect = expectedvalue();
            if(array[offset] != expect) {
                fprintf(stderr,"data mismatch: array[%d]=%f expected=%f\n",
                            offset,array[offset],expect);
                errs++;
                if(errs >= MAXERRS)
                    break;
            }
            odom_next();
       }
   }

   if(errs == 0)
        fprintf(stderr,"no data errors\n");
   return (errs == 0);
}

static void
showparameters(void)
{
    int i;
    fprintf(stderr,"test: nparams=%ld: params=",(unsigned long)nparams);
    for(i=0;i<nparams;i++) {
        fprintf(stderr," %u",params[i]);
    }
    fprintf(stderr,"\n");
    for(i=0;i<actualdims;i++)
        fprintf(stderr,"%s%ld",(i==0?" chunks=":","),(unsigned long)chunks[i]);
    fprintf(stderr,"\n");
    fflush(stderr);
}

static void
insert(int index, void* src, size_t size)
{
    void* dst = &baseline[index];
    memcpy(dst,src,size);
}

static void
buildbaseline(unsigned int testcasenumber)
{
    unsigned int val4;
    unsigned long long val8;
    float float4;
    double float8;

    baseline[0] = testcasenumber;
    switch (testcasenumber) {
    case 1:
        val4 = ((unsigned int)-17) & 0xff;
        insert(1,&val4,sizeof(val4)); /* 1 signed int*/
	val4 = (unsigned int)23;
        insert(2,&val4,sizeof(val4)); /* 2 unsigned int*/
        val4 = ((unsigned int)-25) & 0xffff;
        insert(3,&val4,sizeof(val4)); /* 3 signed int*/
	val4 = (unsigned int)27;
        insert(4,&val4,sizeof(val4)); /* 4 unsigned int*/
	val4 = (unsigned int)77;
        insert(5,&val4,sizeof(val4)); /* 5 signed int*/
	val4 = (unsigned int)93;
        insert(6,&val4,sizeof(val4)); /* 6 unsigned int*/
	float4 = 789.0f;
        insert(7,&float4,sizeof(float4)); /* 7 float */
	float8 = DBLVAL;
        insert(8,&float8,sizeof(float8)); /* 8 double */
	val8 = -9223372036854775807L;
        insert(10,&val8,sizeof(val8)); /* 10 signed long long */
	val8 = 18446744073709551615UL;
        insert(12,&val8,sizeof(val8)); /* 12 unsigned long long */
	break;
    default:
	fprintf(stderr,"Unknown testcase number: %d\n",testcasenumber);
	abort();
    }
}

static int
test_test1(void)
{
    int ok = 1;

    reset();

    buildbaseline(1);

    fprintf(stderr,"test1: compression.\n");
    create();
    setchunking();
    setvarfilter();
    showparameters();
    CHECK(nc_enddef(ncid));

    /* Fill in the array */
    fill();
    /* write array */
    CHECK(nc_put_var(ncid,varid,expected));
    CHECK(nc_close(ncid));

    fprintf(stderr,"test1: decompression.\n");
    reset();
    openfile();
    CHECK(nc_get_var_float(ncid, varid, array));
    ok = compare();
    CHECK(nc_close(ncid));
    return ok;
}

/**************************************************/
/* Utilities */

static void
reset()
{
    memset(array,0,sizeof(float)*actualproduct);
}

static void
odom_reset(void)
{
    memset(odom,0,sizeof(odom));
}

static int
odom_more(void)
{
    return (odom[0] < dims[0]);
}

static int
odom_next(void)
{
    int i; /* do not make unsigned */
    for(i=actualdims-1;i>=0;i--) {
        odom[i] += 1;
        if(odom[i] < dims[i]) break;
        if(i == 0) return 0; /* leave the 0th entry if it overflows*/
        odom[i] = 0; /* reset this position*/
    }
    return 1;
}

static int
odom_offset(void)
{
    int i;
    int offset = 0;
    for(i=0;i<actualdims;i++) {
        offset *= dims[i];
        offset += odom[i];
    }
    return offset;
}

static float
expectedvalue(void)
{
    int i;
    float offset = 0;

    for(i=0;i<actualdims;i++) {
        offset *= dims[i];
        offset += odom[i];
    }
    return offset;
}

static void
init(int argc, char** argv)
{
    int i;
    /* Setup various variables */
    totalproduct = 1;
    actualproduct = 1;
    chunkproduct = 1;
    for(i=0;i<MAXDIMS;i++) {
        dims[i] = dimsize;
        chunks[i] = (pattern[i] == 1 ? 1 : chunksize);
        totalproduct *= dims[i];
        if(i < actualdims) {
            actualproduct *= dims[i];
            chunkproduct *= chunks[i];
        }
    }
    /* Allocate max size */
    array = (float*)calloc(1,sizeof(float)*actualproduct);
    expected = (float*)calloc(1,sizeof(float)*actualproduct);
}

/**************************************************/
int
main(int argc, char **argv)
{
#ifdef DEBUG
    H5Eprint(stderr);
    nc_set_log_level(1);
#endif
    init(argc,argv);
    if(!test_test1()) ERRR;
    exit(nerrs > 0?1:0);
}
