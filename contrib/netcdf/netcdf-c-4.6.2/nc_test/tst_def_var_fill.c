/* This is part of the netCDF package.
 * Copyright 2005 University Corporation for Atmospheric Research/Unidata
 * See COPYRIGHT file for conditions of use.
 *
 * Test per-variable fill mode for classic file formats.
 *
 * Author: Wei-keng Liao.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <netcdf.h>


#define CHECK_ERR { \
    if (err != NC_NOERR) { \
        nerrs++; \
        printf("Error at line %d in %s: (%s)\n", \
        __LINE__,__FILE__,nc_strerror(err)); \
    } \
}

#define EXP_ERR(exp) { \
    if (err != exp) { \
        nerrs++; \
        printf("Error at line %d in %s: expecting %s but got %s\n", \
        __LINE__,__FILE__,#exp, nc_strerror(err)); \
    } \
}

#define NY 8
#define NX 5

int main(int argc, char** argv) {
    char filename[256];
    int i, j, k, err, nerrs=0, ncid, varid[2], dimid[2], *buf;
    size_t start[2], count[2];
    int formats[5]={NC_FORMAT_CLASSIC, NC_FORMAT_64BIT_OFFSET, NC_FORMAT_CDF5,
                    NC_FORMAT_NETCDF4, NC_FORMAT_NETCDF4_CLASSIC};

    if (argc > 2) {
        printf("Usage: %s [filename]\n",argv[0]);
        return 1;
    }
    if (argc == 2) snprintf(filename, 256, "%s", argv[1]);
    else           strcpy(filename, "tst_def_var_fill.nc");

    char *cmd_str = (char*)malloc(strlen(argv[0]) + 256);
    sprintf(cmd_str, "*** TESTING C   %s for def_var_fill ", argv[0]);
    printf("%-66s ------ ", cmd_str); fflush(stdout);
    free(cmd_str);

    buf = (int*) malloc(NY*NX * sizeof(int));

    for (k=0; k<5; k++) {
#ifndef ENABLE_CDF5
        if (formats[k] == NC_FORMAT_CDF5) continue;
#endif
#ifndef USE_NETCDF4
        if (formats[k] == NC_FORMAT_NETCDF4 ||
            formats[k] == NC_FORMAT_NETCDF4_CLASSIC)
            continue;
#endif
        nc_set_default_format(formats[k], NULL);

        /* create a new file for writing ------------------------------------*/
        err = nc_create(filename, NC_CLOBBER, &ncid); CHECK_ERR

        /* define dimension */
        err = nc_def_dim(ncid, "Y", NY, &dimid[0]); CHECK_ERR
        err = nc_def_dim(ncid, "X", NX, &dimid[1]); CHECK_ERR

        /* define variables */
        err = nc_def_var(ncid, "var_nofill", NC_INT, 2, dimid, &varid[0]); CHECK_ERR
        err = nc_def_var(ncid, "var_fill",   NC_INT, 2, dimid, &varid[1]); CHECK_ERR

        /* set fill mode for variables */
        err = nc_def_var_fill(ncid, NC_GLOBAL, 0, NULL); EXP_ERR(NC_EGLOBAL)
        err = nc_def_var_fill(ncid, varid[0], 1, NULL); CHECK_ERR
        err = nc_def_var_fill(ncid, varid[1], 0, NULL); CHECK_ERR

        err = nc_enddef(ncid); CHECK_ERR

        /* write a subarray to both variables */
        for (i=0; i<NY*NX; i++) buf[i] = 5;
        start[0] = 0;
        start[1] = 2;
        count[0] = NY;
        count[1] = 2;
        err = nc_put_vara_int(ncid, varid[0], start, count, buf); CHECK_ERR
        err = nc_put_vara_int(ncid, varid[1], start, count, buf); CHECK_ERR
        err = nc_close(ncid); CHECK_ERR

        /* Now, reopen the file and read variables back */
        err = nc_open(filename, NC_WRITE, &ncid); CHECK_ERR

        /* get variable IDs */
        err = nc_inq_varid(ncid, "var_nofill", &varid[0]); CHECK_ERR
        err = nc_inq_varid(ncid, "var_fill",   &varid[1]); CHECK_ERR

        /* read variable "var_nofill" and check contents */
        for (i=0; i<NY*NX; i++) buf[i] = -1;
        err = nc_get_var_int(ncid, varid[0], buf); CHECK_ERR
        for (i=0; i<NY; i++) {
            for (j=0; j<NX; j++) {
                if (2 <= j && j < 4) {
                    if (buf[i*NX+j] != 5) {
                        printf("Error at line %d in %s: expect get buf[%d]=%d but got %d\n",
                               __LINE__,__FILE__,i*NX+j, 5, buf[i*NX+j]);
                        nerrs++;
                    }
                }
                else if (buf[i*NX+j] == NC_FILL_INT) {
                    printf("Warning at line %d in %s: get buf[%d] same as NC_FILL_INT\n",
                           __LINE__,__FILE__,i*NX+j);
                }
            }
        }

        /* read variable "var_fill" and check contents */
        for (i=0; i<NY*NX; i++) buf[i] = -1;
        err = nc_get_var_int(ncid, varid[1], buf); CHECK_ERR
        for (i=0; i<NY; i++) {
            for (j=0; j<NX; j++) {
                int expect = NC_FILL_INT;
                if (2 <= j && j< 4) expect = 5;

                if (buf[i*NX+j] != expect) {
                    printf("Error at line %d in %s: expect get buf[%d]=%d but got %d\n",
                           __LINE__,__FILE__,i*NX+j, expect, buf[i*NX+j]);
                    nerrs++;
                }
            }
        }

        err = nc_close(ncid); CHECK_ERR
    }
    free(buf);

    if (nerrs) printf("fail with %d mismatches\n",nerrs);
    else       printf("pass\n");

    return (nerrs > 0);
}

