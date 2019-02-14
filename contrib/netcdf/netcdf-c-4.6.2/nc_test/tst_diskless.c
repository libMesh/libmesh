/** \file \internal
Basic diskless API tests.

Copyright 2011, UCAR/Unidata. See COPYRIGHT file for copying and
redistribution conditions.
*/

#define DDBG

#include "config.h"
#include <stdio.h>
#include <stdlib.h>
#include "netcdf.h"
#include "nc_tests.h"
#include "err_macros.h"

/*
netcdf tst_diskless {
variables:
   int resistor_value;
   float capacitor_value;
   short number_of_555_timer_chips;
}
*/

#ifndef NC_NETCDF3
#define NC_NETCDF3 0
#endif

#define FLAGS4 (NC_NETCDF4|NC_CLASSIC_MODEL)
#define FLAGS3 (NC_NETCDF3)

#define RESISTOR "resistor_value"
#define CAPACITOR "capacitor_value"
#define NUM555 "number_of_555_timer_chips"

#ifdef DDBG
#undef ERR
void fail(int line) {
    fflush(stdout);
    fprintf(stderr,"\nfail: line=%d\n",line);
    fflush(stderr);
    exit(1);
}
#define ERR fail(__LINE__)
#endif

/* Control flags  */
static int flags, persist, usenetcdf4, mmap, diskless;

char*
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

/* Remove a file; do not care if it does not exist */
static void
removefile(int persist,  char* filename)
{
    if(persist) {
	remove(filename);
    }
}

int
main(int argc, char **argv)
{
    int i;
    char* filename = "tst_diskless.nc";

    /* Set defaults */
    persist = 0;
    usenetcdf4 = 0;
    mmap = 0;
    diskless = 0;

    for(i=1;i<argc;i++) {
	if(strcmp(argv[i],"netcdf4")==0) usenetcdf4=1;
	else if(strcmp(argv[i],"persist")==0) persist=1;
	else if(strcmp(argv[i],"mmap")==0) mmap=1;
	else if(strcmp(argv[i],"diskless")==0) diskless=1;
	/* ignore anything not recognized */
    }

#ifndef USE_NETCDF4
    fprintf(stderr,"netcdf-4 not supported; ignored\n");
    usenetcdf4 = 0;
#endif

    /* Invalid combinations */
    if(mmap && diskless) {fprintf(stderr,"Illegal: mmap+diskless\n"); exit(1);};
    if(mmap && usenetcdf4) {fprintf(stderr,"Illegal: mmap+netcdf4\n"); exit(1);};

    flags = usenetcdf4?FLAGS4:FLAGS3;
    if(persist) flags |= NC_PERSIST;
    if(mmap) flags |= NC_MMAP;
    if(diskless) flags |= NC_DISKLESS;

printf("\n*** Testing the diskless|mmap API.\n");
printf("*** testing diskless file with scalar vars...");
{
    int ncid, varid0, varid1, varid2;
    int ndims_in, nvars_in, natts_in, unlimdimid_in;
    char name_in[NC_MAX_NAME + 1];
    nc_type type_in;
    float float_data = 3.14, float_data_in;
    int int_data = 42, int_data_in;
    short short_data = 2, short_data_in;

    removefile(persist,filename);

    /* Create a netCDF file (which exists in memory). */
    if (nc_create(filename, flags, &ncid)) ERR;

    /* Create some variables. */
    if (nc_def_var(ncid, RESISTOR, NC_INT, 0, NULL, &varid0)) ERR;
    if (nc_def_var(ncid, CAPACITOR, NC_FLOAT, 0, NULL, &varid1)) ERR;
    if (nc_def_var(ncid, NUM555, NC_SHORT, 0, NULL, &varid2)) ERR;
    if (nc_enddef(ncid)) ERR;

    /* Write some data to this file. */
    if (nc_put_vara_int(ncid, varid0, NULL, NULL, &int_data)) ERR;
    if (nc_put_vara_float(ncid, varid1, NULL, NULL, &float_data)) ERR;
    if (nc_put_vara_short(ncid, varid2, NULL, NULL, &short_data)) ERR;

    /* Now check the phony file. */
    if (nc_inq(ncid, &ndims_in, &nvars_in, &natts_in, &unlimdimid_in)) ERR;
    if (ndims_in != 0 || nvars_in != 3 || natts_in != 0 || unlimdimid_in != -1) ERR;

    /* Check variables. */
    if (nc_inq_var(ncid, varid0, name_in, &type_in, &ndims_in, NULL, &natts_in)) ERR;
    if (strcmp(name_in, RESISTOR) || type_in != NC_INT || ndims_in != 0 ||
    natts_in != 0) ERR;
    if (nc_inq_var(ncid, varid1, name_in, &type_in, &ndims_in, NULL, &natts_in)) ERR;
    if (strcmp(name_in, CAPACITOR) || type_in != NC_FLOAT || ndims_in != 0 ||
    natts_in != 0) ERR;
    if (nc_inq_var(ncid, varid2, name_in, &type_in, &ndims_in, NULL, &natts_in)) ERR;
    if (strcmp(name_in, NUM555) || type_in != NC_SHORT || natts_in != 0) ERR;

    /* Read my absolutely crucial data. */
    if (nc_get_vara_int(ncid, varid0, NULL, NULL, &int_data_in)) ERR;
    if (int_data_in != int_data) ERR;
    if (nc_get_vara_float(ncid, varid1, NULL, NULL, &float_data_in)) ERR;
    if (float_data_in != float_data) ERR;
    if (nc_get_vara_short(ncid, varid2, NULL, NULL, &short_data_in)) ERR;
    if (short_data_in != short_data) ERR;

    /* Close the file. */
    if (nc_close(ncid))
      abort(); /* ERR; */
    }
    SUMMARIZE_ERR;

    if(!usenetcdf4 && persist) {
        int ncid, varid0, varid1, varid2;
        float float_data = 3.14, float_data_in;
        int int_data = 42, int_data_in;
        short short_data = 2, short_data_in;

        printf("*** testing diskless open of previously created file...");

        if (nc_open(filename, flags, &ncid)) ERR;

	/* Read and compare */
        if (nc_inq_varid(ncid, RESISTOR, &varid0)) ERR;
        if (nc_inq_varid(ncid, CAPACITOR, &varid1)) ERR;
        if (nc_inq_varid(ncid, NUM555, &varid2)) ERR;

        if (nc_get_vara_int(ncid, varid0, NULL, NULL, &int_data_in)) ERR;
        if (int_data_in != int_data) ERR;
        if (nc_get_vara_float(ncid, varid1, NULL, NULL, &float_data_in)) ERR;
        if (float_data_in != float_data) ERR;
        if (nc_get_vara_short(ncid, varid2, NULL, NULL, &short_data_in)) ERR;
        if (short_data_in != short_data) ERR;

	nc_close(ncid);
    }
    SUMMARIZE_ERR;

    printf("*** testing creation of simple diskless file...");
    {
    #define NDIMS 2
    #define DIM0_NAME "Fun"
    #define DIM1_NAME "Money"
    #define DIM1_LEN 200
    #define ATT0_NAME "home"
    #define ATT0_TEXT "earthship"
    #define VAR0_NAME "nightlife"
    #define VAR1_NAME "time"
    #define VAR2_NAME "taxi_distance"

    int ncid, dimid[NDIMS], dimid_in[NDIMS], varid0, varid1, varid2;
    int ndims_in, nvars_in, natts_in, unlimdimid_in;
    char name_in[NC_MAX_NAME + 1], att0_in[NC_MAX_NAME + 1];
    nc_type type_in;
    size_t len_in;
    short short_data[DIM1_LEN];
    size_t start[1] = {0};
    size_t count[1] = {DIM1_LEN};
    int i;
    float float_data = 42.22, float_data_in;

    /* This is some really important data that I want to save. */
    for (i = 0; i < DIM1_LEN; i++)
    short_data[i] = i;

    removefile(persist,filename);

    /* Create a netCDF file (which exists only in memory). I am
    * confident that the world-famous netCDF format is the way to
    * store my data! */
    if (nc_create(filename, flags, &ncid)) ERR;

    /* Create some atts. They will help document my data forever. */
    if (nc_put_att_text(ncid, NC_GLOBAL, ATT0_NAME,
			sizeof(ATT0_TEXT), ATT0_TEXT)) ERR;

    /* Create dimensions: money is limited, but fun is not! */
    if (nc_def_dim(ncid, DIM0_NAME, NC_UNLIMITED, &dimid[0])) ERR;
    if (nc_def_dim(ncid, DIM1_NAME, DIM1_LEN, &dimid[1])) ERR;

    /* Create some variables. The data they hold must persist
    * through the ages. */
    if (nc_def_var(ncid, VAR0_NAME, NC_INT, NDIMS, dimid, &varid0)) ERR;
    if (nc_def_var(ncid, VAR1_NAME, NC_FLOAT, 0, NULL, &varid1)) ERR;
    if (nc_def_var(ncid, VAR2_NAME, NC_SHORT, 1, &dimid[1], &varid2)) ERR;
    if (nc_enddef(ncid)) ERR;

    /* Write some data to this file. I'm glad I'm saving this
    * important data in such a safe format! */
    if (nc_put_vara_float(ncid, varid1, NULL, NULL, &float_data)) ERR;
    if (nc_put_vara_short(ncid, varid2, start, count, short_data)) ERR;

    /* Now check the phony file. Is my data safe? */
    if (nc_inq(ncid, &ndims_in, &nvars_in, &natts_in, &unlimdimid_in)) ERR;
    if (ndims_in != 2 || nvars_in != 3 || natts_in != 1 || unlimdimid_in != 0) ERR;

    /* Check attributes - they will be needed by future generations
    * of scientists to understand my data. */
    if (nc_get_att_text(ncid, NC_GLOBAL, ATT0_NAME, att0_in)) ERR;
    att0_in[sizeof(ATT0_TEXT)] = '\0';
    if (strcmp(att0_in, ATT0_TEXT)) ERR;

    /* Check dimensions. */
    if (nc_inq_dim(ncid, dimid[0], name_in, &len_in)) ERR;
    if (strcmp(name_in, DIM0_NAME) || len_in != 0) ERR;
    if (nc_inq_dim(ncid, dimid[1], name_in, &len_in)) ERR;
    if (strcmp(name_in, DIM1_NAME) || len_in != DIM1_LEN) ERR;

    /* Check variables. */
    if (nc_inq_var(ncid, varid0, name_in, &type_in, &ndims_in, dimid_in, &natts_in)) ERR;
    if (strcmp(name_in, VAR0_NAME) || type_in != NC_INT || ndims_in != NDIMS ||
    dimid_in[0] != 0 || dimid_in[1] != 1 || natts_in != 0) ERR;
    if (nc_inq_var(ncid, varid1, name_in, &type_in, &ndims_in, dimid_in, &natts_in)) ERR;
    if (strcmp(name_in, VAR1_NAME) || type_in != NC_FLOAT || ndims_in != 0 ||
    natts_in != 0) ERR;
    if (nc_inq_var(ncid, varid2, name_in, &type_in, &ndims_in, dimid_in, &natts_in)) ERR;
    if (strcmp(name_in, VAR2_NAME) || type_in != NC_SHORT || ndims_in != 1 ||
    dimid_in[0] != 1 || natts_in != 0) ERR;

    /* Read my absolutely crucial data. */
    if (nc_get_vara_float(ncid, varid1, NULL, NULL, &float_data_in)) ERR;
    if (float_data_in != float_data) ERR;

    /* Close the file, losing all information. Hey! What kind of
    * storage format is this, anyway? */
    if (nc_close(ncid))
      abort(); /* ERR; */
    }
    SUMMARIZE_ERR;
    printf("*** testing diskless file with scalar vars and type conversion...");
    {
    #define DUNE "dune"
    #define STAR_TREK "capacitor_value"
    #define STAR_WARS "number_of_555_timer_chips"

    int ncid, varid0, varid1, varid2;
    int ndims_in, nvars_in, natts_in, unlimdimid_in;
    char name_in[NC_MAX_NAME + 1];
    nc_type type_in;
    float float_data = 3.14, float_data_in;
    int int_data = 42, int_data_in;
    short short_data = 2, short_data_in;

    removefile(persist,filename);

    /* Create a netCDF file (which exists only in memory). */
    if (nc_create(filename, flags, &ncid)) ERR;

    /* Create some variables. */
    if (nc_def_var(ncid, DUNE, NC_INT, 0, NULL, &varid0)) ERR;
    if (nc_def_var(ncid, STAR_TREK, NC_FLOAT, 0, NULL, &varid1)) ERR;
    if (nc_def_var(ncid, STAR_WARS, NC_SHORT, 0, NULL, &varid2)) ERR;
    if (nc_enddef(ncid)) ERR;

    /* Write some data to this file. */
    if (nc_put_vara_int(ncid, varid0, NULL, NULL, &int_data)) ERR;
    if (nc_put_vara_float(ncid, varid1, NULL, NULL, &float_data)) ERR;
    if (nc_put_vara_short(ncid, varid2, NULL, NULL, &short_data)) ERR;

    /* Now check the phony file. */
    if (nc_inq(ncid, &ndims_in, &nvars_in, &natts_in, &unlimdimid_in)) ERR;
    if (ndims_in != 0 || nvars_in != 3 || natts_in != 0 || unlimdimid_in != -1) ERR;

    /* Check variables. */
    if (nc_inq_var(ncid, varid0, name_in, &type_in, &ndims_in, NULL, &natts_in)) ERR;
    if (strcmp(name_in, DUNE) || type_in != NC_INT || ndims_in != 0 ||
    natts_in != 0) ERR;
    if (nc_inq_var(ncid, varid1, name_in, &type_in, &ndims_in, NULL, &natts_in)) ERR;
    if (strcmp(name_in, STAR_TREK) || type_in != NC_FLOAT || ndims_in != 0 ||

    natts_in != 0) ERR;
    if (nc_inq_var(ncid, varid2, name_in, &type_in, &ndims_in, NULL, &natts_in)) ERR;
    if (strcmp(name_in, STAR_WARS) || type_in != NC_SHORT || natts_in != 0) ERR;

    /* Read my absolutely crucial data. */
    if (nc_get_vara_int(ncid, varid0, NULL, NULL, &int_data_in)) ERR;
    if (int_data_in != int_data) ERR;
    if (nc_get_vara_float(ncid, varid1, NULL, NULL, &float_data_in)) ERR;
    if (float_data_in != float_data) ERR;
    if (nc_get_vara_short(ncid, varid2, NULL, NULL, &short_data_in)) ERR;
    if (short_data_in != short_data) ERR;

    /* Close the file. */
    if (nc_close(ncid))
      abort(); /* ERR; */
    }
    SUMMARIZE_ERR;
    FINAL_RESULTS;

    /* Unnecessary exit(0), FINAL_RESULTS returns. */
    /* exit(0); */
}
