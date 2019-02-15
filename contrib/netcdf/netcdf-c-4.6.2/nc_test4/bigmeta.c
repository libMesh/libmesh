/*
Create a netcdf-4 file with horrendously large metadata.
*/

#include "config.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#include <netcdf.h>
#ifdef HAVE_GETOPT_H
#include <getopt.h>
#endif

#undef DEBUG

/*
Total number of groups is NGROUPS^TREEDEPTH
TREEDEPTH =  Depth of the group tree
NGROUPS =  Number subgroups per group
NGROUPATTRS =  Number group attributes per group
NDIMS =  Number dimensions per group
NTYPES =  Number user types per group
NVARS =  Number of variables per group
VARRANK =  Rank of variables: must be <= NDIMS
NVARATTRS =  Number attributes per variable
*/

/* Define the Defaults */
#if 1
#define TREEDEPTH 6
#define NGROUPS 2
#define NGROUPATTRS 100
#define NDIMS 100
#define NTYPES 10
#define NVARS 100
#define VARRANK 2
#define NVARATTRS 500
#endif

/* Alternate Defaults for testing*/
#if 0
#define TREEDEPTH 2
#define NGROUPS 2
#define NGROUPATTRS 2
#define NDIMS 2
#define NTYPES 2
#define NVARS 2
#define VARRANK 2
#define NVARATTRS NGROUPATTRS
#endif

#define FILE "bigmeta.nc"

#define CHECK(expr) assert((expr) == NC_NOERR)

static int treedepth = TREEDEPTH;
static int ngroups = NGROUPS;
static int ngroupattrs = NGROUPATTRS;
static int ndims = NDIMS;
static int ntypes = NTYPES;
static int nvars = NVARS;
static int varrank = VARRANK;
static int nvarattrs = NVARATTRS;

/* Define the getopt tags */
#define OPT_UNKNOWN 0
#define OPT_TREEDEPTH 1
#define OPT_NGROUPS 2
#define OPT_NGROUPATTRS 3
#define OPT_NDIMS 4
#define OPT_NTYPES 5
#define OPT_NVARS 6
#define OPT_VARRANK 7
#define OPT_NVARATTRS 8

static struct option options[] = {
{"treedepth", 1, NULL, OPT_TREEDEPTH},
{"ngroups", 1, NULL, OPT_NGROUPS},
{"ngroupattrs", 1, NULL, OPT_NGROUPATTRS},
{"ndims", 1, NULL, OPT_NDIMS},
{"ntypes", 1, NULL, OPT_NTYPES},
{"nvars", 1, NULL, OPT_NVARS},
{"varrank", 1, NULL, OPT_VARRANK},
{"nvarattrs", 1, NULL, OPT_NVARATTRS},
{NULL, 0, NULL, 0}
};

/**************************************************/

static void
reportparameters(void)
{
    fprintf(stderr,"--treedepth=%d\n",treedepth);
    fprintf(stderr,"--ngroups=%d\n",ngroups);
    fprintf(stderr,"--ngroupattrs=%d\n",ngroupattrs);
    fprintf(stderr,"--ndims=%d\n",ndims);
    fprintf(stderr,"--ntypes=%d\n",ntypes);
    fprintf(stderr,"--nvars=%d\n",nvars);
    fprintf(stderr,"--varrank=%d\n",varrank);
    fprintf(stderr,"--nvarattrs=%d\n",nvarattrs);
    fflush(stderr);
}

/* Build a compound type with two fields */
static void
buildcmpdtype(int grpid, int typindex)
{
    char name[NC_MAX_NAME+1];
    int typid;
    snprintf(name,NC_MAX_NAME,"cmpd%d",typindex);
    CHECK(nc_def_compound(grpid, 2*sizeof(int), name, &typid));
    CHECK(nc_insert_compound(grpid, typid, "f1", 0, NC_INT));
    CHECK(nc_insert_compound(grpid, typid, "f2", sizeof(int), NC_INT));
}

/* Build an enum type with two members */
static void
buildenumtype(int grpid, int typindex)
{
    char name[NC_MAX_NAME+1];
    int typid;
    int val0 = 17;
    int val1 = 37;
    snprintf(name,NC_MAX_NAME,"enum%d",typindex);
    CHECK(nc_def_enum(grpid, NC_INT, name, &typid));
    CHECK(nc_insert_enum(grpid, typid, "m0", &val0));
    CHECK(nc_insert_enum(grpid, typid, "m1", &val1));
}

/* Build attributes for either group or var */
static void
buildatts(int grpid, int varid)
{
    char name[NC_MAX_NAME+1];
    int i, count;

    count = (varid == NC_GLOBAL? ngroupattrs : nvarattrs);

    for(i=0;i<count;i++) {
        snprintf(name,NC_MAX_NAME,"a%d",i);
        CHECK(nc_put_att_int(grpid, varid, name, NC_INT, 1, &i));
    }
}

static void
buildgroup(int parent, int grpindex, int depth)
{
    char name[NC_MAX_NAME+1];
    int i, grpid, varid;
    int dimids[NDIMS];

    if(depth == 0) return;

    snprintf(name,NC_MAX_NAME,"g%d",grpindex);
    CHECK(nc_def_grp(parent, name, &grpid));
 
    /* Add dimensions and capture ids */
    for(i=0;i<ndims;i++) {
	snprintf(name,NC_MAX_NAME,"d%d",i);
        CHECK(nc_def_dim(grpid, name, i, &dimids[i]));
    }

    /* Add types: even index => compound, odd => enum */
    for(i=0;i<ntypes;i++) {
	if((i % 2) == 0)
	    buildcmpdtype(grpid,i);
	else
	    buildenumtype(grpid,i);
    }

    /* Add variables (plus their attributes */
    for(i=0;i<nvars;i++) {
	snprintf(name,NC_MAX_NAME,"v%d",i);
        CHECK(nc_def_var(grpid, name, NC_INT, VARRANK, dimids, &varid));
	/* Make variable be nofill */
	CHECK(nc_def_var_fill(grpid,varid,1,NULL));
	buildatts(grpid,varid);
    }

    /* Add group attributes */
    buildatts(grpid,NC_GLOBAL);

    /* Build subgroups */
    for(i=0;i<ngroups;i++) {
	buildgroup(grpid,i,depth-1);
    }

}

int
main(int argc, char **argv)
{
    int i, ncid;
    time_t starttime, endtime;
    long long delta;
    int tag;

    if(argc > 1) {
	while ((tag = getopt_long_only(argc, argv, "", options, NULL)) >= 0) {
#ifdef DEBUG
fprintf(stderr,"arg=%s value=%s\n",argv[optind-1],optarg);
#endif
	    switch (tag) {
	    case OPT_TREEDEPTH:
		treedepth = atoi(optarg);
		break;
	    case OPT_NGROUPS:
		ngroups = atoi(optarg);
		break;
	    case OPT_NGROUPATTRS:
		ngroupattrs = atoi(optarg);
		break;
	    case OPT_NDIMS:
		ndims = atoi(optarg);
		break;
	    case OPT_NTYPES:
		ntypes = atoi(optarg);
		break;
	    case OPT_NVARS:
		nvars = atoi(optarg);
		break;
	    case OPT_VARRANK:
		varrank = atoi(optarg);
		break;
	    case OPT_NVARATTRS:
		nvarattrs = atoi(optarg);
		break;
	    case ':':
		fprintf(stderr,"missing argument\n");
		exit(1);
	    case '?':
	    default:
		fprintf(stderr,"unknown option\n");
		exit(1);
	    }
	}
    }

    reportparameters();

    starttime = 0;
    endtime = 0;

    time(&starttime);

    CHECK(nc_create(FILE, NC_NETCDF4, &ncid));
    /* Build subgroups */
    for(i=0;i<NGROUPS;i++) {
	buildgroup(ncid,i,TREEDEPTH - 1);
    }

    CHECK(nc_close(ncid));

    time(&endtime);

    /* Compute the delta 1 second resolution is fine for this */
    delta = (long long)(endtime - starttime);
    printf("create delta=%lld\n",delta);
    return 0;

    return 0;
}
