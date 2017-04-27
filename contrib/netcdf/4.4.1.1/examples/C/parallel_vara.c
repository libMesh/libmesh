/*********************************************************************
 *
 *  Copyright (C) 2014, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *********************************************************************/

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * This example shows how to use nc_put_vara_int() to write a 2D 4-byte integer
 * array in parallel and read it back using the same array partitioning pattern.
 * It first defines a netCDF variable of size global_nx * global_ny where
 *    global_ny == NY and
 *    global_nx == (NX * number of MPI processes).
 * The data partitioning pattern is a column-wise partitioning across all
 * proceses. Each process writes a subarray of size ny * nx.
 *
 *    To compile:
 *        mpicc -O2 parallel_vara.c -o parallel_vara -lnetcdf -lpnetcdf
 *
 * Example commands for MPI run and outputs from running ncdump on the
 * NC file produced by this example program:
 *
 *    % mpiexec -n 4 ./parallel_vara /pvfs2/wkliao/testfile.nc
 *
 *    % ncdump /pvfs2/wkliao/testfile.nc
 *    netcdf testfile {
 *    dimensions:
 *            y = 10 ;
 *            x = 16 ;
 *    variables:
 *            int var(y, x) ;
 *                var:str_att_name = "example attribute of type text." ;
 *                var:float_att_name = 0.f, 1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 7.f ;
 *    // global attributes:
 *                :history = "Wed Apr 30 11:18:58 2014\n",
 *       "" ;
 *    data:
 *
 *     var =
 *         0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3,
 *         0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3,
 *         0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3,
 *         0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3,
 *         0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3,
 *         0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3,
 *         0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3,
 *         0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3,
 *         0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3,
 *         0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3 ;
 *    }
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* strcpy() */
#include <unistd.h> /* getopt() */
#include <time.h>   /* time() localtime(), asctime() */
#include <assert.h>
#include <mpi.h>
#include <netcdf.h>
#include <netcdf_par.h>

#define NY 10
#define NX 4

#define FATAL_ERR {if(err!=NC_NOERR) {printf("Error at line=%d: %s Aborting ...\n", __LINE__, nc_strerror(err)); goto fn_exit;}}
#define ERR {if(err!=NC_NOERR)printf("Error at line=%d: %s\n", __LINE__, nc_strerror(err));}

static void
usage(char *argv0)
{
    char *help =
    "Usage: %s [-h] | [-q] [file_name]\n"
    "       [-h] Print help\n"
    "       [-q] Quiet mode (reports when fail)\n"
    "       [filename] output netCDF file name\n";
    fprintf(stderr, help, argv0);
}

int main(int argc, char** argv)
{
    extern int optind;
    char filename[128];
    int i, j, rank, nprocs, verbose=1, err;
    int ncid, cmode, omode, varid, dimid[2], buf[NY][NX];
    char str_att[128];
    float float_att[100];
    size_t global_ny, global_nx, start[2], count[2];

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    /* get command-line arguments */
    while ((i = getopt(argc, argv, "hq")) != EOF)
        switch(i) {
            case 'q': verbose = 0;
                      break;
            case 'h':
            default:  if (rank==0) usage(argv[0]);
                      MPI_Finalize();
                      return 0;
        }
    argc -= optind;
    argv += optind;
    if (argc == 1) strcpy(filename, argv[0]); /* optional argument */
    else strcpy(filename, "testfile.nc");

    MPI_Bcast(filename, 128, MPI_CHAR, 0, MPI_COMM_WORLD);

    /* create a new file for writing ----------------------------------------*/
    cmode = NC_CLOBBER | NC_PNETCDF;
    err = nc_create_par(filename, cmode, MPI_COMM_WORLD, MPI_INFO_NULL, &ncid); FATAL_ERR

    /* the global array is NY * (NX * nprocs) */
    global_ny = NY;
    global_nx = NX * nprocs;

    for (i=0; i<NY; i++)
        for (j=0; j<NX; j++)
             buf[i][j] = rank;

    /* add a global attribute: a time stamp at rank 0 */
    time_t ltime = time(NULL); /* get the current calendar time */
    asctime_r(localtime(&ltime), str_att);

    /* make sure the time string are consistent among all processes */
    MPI_Bcast(str_att, strlen(str_att), MPI_CHAR, 0, MPI_COMM_WORLD);

    err = nc_put_att_text(ncid, NC_GLOBAL, "history", strlen(str_att),
                          &str_att[0]); ERR

    /* define dimensions x and y */
    err = nc_def_dim(ncid, "Y", global_ny, &dimid[0]); ERR
    err = nc_def_dim(ncid, "X", global_nx, &dimid[1]); ERR

    /* define a 2D variable of integer type */
    err = nc_def_var(ncid, "var", NC_INT, 2, dimid, &varid); ERR

    /* add attributes to the variable */
    strcpy(str_att, "example attribute of type text.");
    err = nc_put_att_text(ncid, varid, "str_att_name", strlen(str_att),
                          &str_att[0]); ERR

    for (i=0; i<8; i++) float_att[i] = i;
    err = nc_put_att_float(ncid, varid, "float_att_name", NC_FLOAT, 8,
                           &float_att[0]); ERR

    /* do not forget to exit define mode */
    err = nc_enddef(ncid); ERR

    /* set to use MPI/PnetCDF collective I/O */
    err = nc_var_par_access(ncid, varid, NC_COLLECTIVE); ERR

    /* now we are in data mode */
    start[0] = 0;
    start[1] = NX * rank;
    count[0] = NY;
    count[1] = NX;

    err = nc_put_vara_int(ncid, varid, start, count, &buf[0][0]); ERR

    err = nc_close(ncid); ERR

    omode = NC_PNETCDF | NC_NOWRITE;
    err = nc_open_par(filename, omode, MPI_COMM_WORLD, MPI_INFO_NULL, &ncid); FATAL_ERR

    /* inquire dimension IDs and lengths */
    err = nc_inq_dimid(ncid, "Y", &dimid[0]); ERR
    err = nc_inq_dimid(ncid, "X", &dimid[1]); ERR

    err = nc_inq_dimlen(ncid, dimid[0], &global_ny); ERR
    err = nc_inq_dimlen(ncid, dimid[1], &global_nx); ERR

    /* obtain variable ID */
    err = nc_inq_varid(ncid, "var", &varid); ERR

    /* set to use MPI/PnetCDF collective I/O */
    err = nc_var_par_access(ncid, varid, NC_COLLECTIVE); ERR

    /* each process reads its subarray from the file */
    err = nc_get_vara_int(ncid, varid, start, count, &buf[0][0]); ERR

    /* close the file */
    err = nc_close(ncid); ERR

fn_exit:
    MPI_Finalize();
    return 0;
}

