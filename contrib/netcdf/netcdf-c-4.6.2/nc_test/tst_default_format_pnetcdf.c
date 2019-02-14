#include "config.h"
#include <stdio.h>
#include <stdlib.h>
#include <netcdf.h>
#include <netcdf_par.h>

#define ERR { \
    if (err != exp_err) { \
        printf("Error at %s line %d: %s\n", __FILE__, __LINE__, \
               nc_strerror(err)); \
        nerrs++; \
    } \
}

static int default_format;

static int
create_check_pnetcdf(char *fname, int cmode, int exp_format)
{
    int nerrs=0, err, exp_err=NC_NOERR, ncid, format;
    char *exp_str;

#ifndef USE_PNETCDF
    exp_err = NC_ENOTBUILT;
#endif

    switch (exp_format) {
        case NC_FORMAT_CLASSIC:      exp_str="NC_FORMAT_CLASSIC";      break;
        case NC_FORMAT_64BIT_OFFSET: exp_str="NC_FORMAT_64BIT_OFFSET"; break;
        case NC_FORMAT_64BIT_DATA:   exp_str="NC_FORMAT_64BIT_DATA";   break;
        case NC_FORMAT_NETCDF4:      exp_str="NC_FORMAT_NETCDF4";      break;
        case NC_FORMAT_NETCDF4_CLASSIC:
                                     exp_str="NC_FORMAT_NETCDF4_CLASSIC";break;
        default: break;
    }

#ifndef ENABLE_CDF5
    if (cmode & NC_64BIT_DATA) exp_err = NC_ENOTBUILT;
#endif

    /* create a file */
    cmode |= NC_CLOBBER;
    err = nc_create_par(fname, cmode, MPI_COMM_WORLD, MPI_INFO_NULL, &ncid); ERR
    if (exp_err == NC_ENOTBUILT) return 0;

    err = nc_close(ncid); ERR

    /* open the file and check its format */
    err = nc_open(fname, NC_NOWRITE, &ncid); ERR
    err = nc_inq_format(ncid, &format); ERR
    if (format != exp_format) {
        char *f_str="", *d_str="";
        switch (format) {
            case NC_FORMAT_CLASSIC:      f_str = "NC_FORMAT_CLASSIC";
                                         break;
            case NC_FORMAT_64BIT_OFFSET: f_str = "NC_FORMAT_64BIT_OFFSET";
                                         break;
            case NC_FORMAT_64BIT_DATA:   f_str = "NC_FORMAT_64BIT_DATA";
                                         break;
            case NC_FORMAT_NETCDF4:      f_str = "NC_FORMAT_NETCDF4";
                                         break;
            case NC_FORMAT_NETCDF4_CLASSIC: f_str = "NC_FORMAT_NETCDF4_CLASSIC";
                                         break;
            default: break;
        }
        switch (default_format) {
            case NC_FORMAT_CLASSIC:      d_str = "NC_FORMAT_CLASSIC";
                                         break;
            case NC_FORMAT_64BIT_OFFSET: d_str = "NC_FORMAT_64BIT_OFFSET";
                                         break;
            case NC_FORMAT_64BIT_DATA:   d_str = "NC_FORMAT_64BIT_DATA";
                                         break;
            case NC_FORMAT_NETCDF4:      d_str = "NC_FORMAT_NETCDF4";
                                         break;
            case NC_FORMAT_NETCDF4_CLASSIC: d_str = "NC_FORMAT_NETCDF4_CLASSIC";
                                         break;
            default: break;
        }

        printf("Error at %s line %d: default is %s and expect %s but got %s\n",
               __FILE__, __LINE__, d_str, exp_str, f_str);
        nerrs++;
    }
    err = nc_close(ncid); ERR
    return nerrs;
}

int main(int argc, char *argv[])
{
    char *fname="tst_default_format.nc";
    int err, exp_err=NC_NOERR, nerrs=0, ncid, cmode;

    MPI_Init(&argc, &argv);

    if (argc == 2) fname = argv[1];

    default_format = NC_FORMAT_CLASSIC;

    /* check illegal cmode */
    cmode = NC_64BIT_OFFSET | NC_64BIT_DATA;
    err = nc_create_par(fname, cmode, MPI_COMM_WORLD, MPI_INFO_NULL, &ncid);
    if (err != NC_EINVAL) {
        printf("Error at %s line %d: expect NC_EINVAL but got %d\n",
               __FILE__, __LINE__, err);
        nerrs++;
    }
#ifdef USE_NETCDF4
    /* check illegal cmode */
    cmode = NC_NETCDF4 | NC_64BIT_OFFSET;
    err = nc_create_par(fname, cmode, MPI_COMM_WORLD, MPI_INFO_NULL, &ncid);
    if (err != NC_EINVAL) {
        printf("Error at %s line %d: expect NC_EINVAL but got %d\n",
               __FILE__, __LINE__, err);
        nerrs++;
    }
#endif

    /* create a file in CDF1 format */
    cmode = 0;
    nerrs += create_check_pnetcdf(fname, cmode, NC_FORMAT_CLASSIC);

    /* create a file in CDF2 format */
    cmode = NC_64BIT_OFFSET;
    nerrs += create_check_pnetcdf(fname, cmode, NC_FORMAT_64BIT_OFFSET);

#ifdef ENABLE_CDF5
    /* create a file in CDF5 format */
    cmode = NC_64BIT_DATA;
    nerrs += create_check_pnetcdf(fname, cmode, NC_FORMAT_64BIT_DATA);
#endif

    /* set default file format to NC_FORMAT_64BIT_OFFSET ------------------*/
    default_format = NC_FORMAT_64BIT_OFFSET;
    err = nc_set_default_format(default_format, NULL); ERR

    /* create a file in default format */
    cmode = 0;
    nerrs += create_check_pnetcdf(fname, cmode, NC_FORMAT_64BIT_OFFSET);

#ifdef ENABLE_CDF5
    /* create a file in CDF5 format (this should ignore default) */
    cmode = NC_64BIT_DATA;
    nerrs += create_check_pnetcdf(fname, cmode, NC_FORMAT_64BIT_DATA);

    /* set default file format to NC_FORMAT_64BIT_DATA --------------------*/
    default_format = NC_FORMAT_64BIT_DATA;
    err = nc_set_default_format(default_format, NULL); ERR

    /* create a file in default format */
    cmode = 0;
    nerrs += create_check_pnetcdf(fname, cmode, NC_FORMAT_64BIT_DATA);
#endif

    /* create a file in CDF2 format (this should ignore default) */
    cmode = NC_64BIT_OFFSET;
    nerrs += create_check_pnetcdf(fname, cmode, NC_FORMAT_64BIT_OFFSET);

    MPI_Finalize();
    return (nerrs > 0);
}
