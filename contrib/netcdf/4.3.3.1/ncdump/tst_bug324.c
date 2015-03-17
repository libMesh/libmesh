#include <nc_tests.h>
#include <stdio.h>
#include <stdlib.h>
#include <netcdf.h>

#define FILENAME "tst_bug324.nc"
#define RANK_LAT 1
#define RANK_H   1
#define LEN_LAT  2
#define LEN_H    2
#define NAME_LAT "lat"
#define NAME_H   "h"

int
main(int argc, char **argv) 
{/* Test bug fix for NCF-324, file that caused nc_close() failure for
  * non-coordinate variable and dimension with the same name */

    int ncid;
    int lat_dim;		/* dimension with associated coordinate variable */
    int h_dim;			/* dimension with no associated coordinate variable */

    size_t lat_len = LEN_LAT;
    size_t h_len = LEN_H;
    int lat_id;
    int h_id;

#   define RANK_LAT 1
#   define RANK_H 1
    int lat_dims[RANK_LAT];
    int h_dims[RANK_H];

   printf("\n*** Testing fix for non-coord var bug.\n");
   printf("*** creating bug test file %s...", FILENAME);

    if (nc_create(FILENAME, NC_CLOBBER|NC_NETCDF4|NC_CLASSIC_MODEL, &ncid)) ERR;
    if (nc_def_dim(ncid, NAME_LAT, lat_len, &lat_dim)) ERR;
    if (nc_def_dim(ncid, NAME_H, h_len, &h_dim)) ERR;
    lat_dims[0] = lat_dim;
    if (nc_def_var(ncid, NAME_LAT, NC_DOUBLE, RANK_LAT, lat_dims, &lat_id)) ERR;
    h_dims[0] = lat_dim;
    if (nc_def_var(ncid, NAME_H, NC_DOUBLE, RANK_H, h_dims, &h_id)) ERR;
    if (nc_enddef (ncid)) ERR;

    {
	double lat_data[LEN_LAT] = {((double)-45), ((double)45)} ;
	size_t lat_startset[1] = {0} ;
	size_t lat_countset[1] = {LEN_LAT};
	if ( nc_put_vara(ncid, lat_id, lat_startset, lat_countset, lat_data) ) ERR;
    }

    {
	double h_data[2] = {((double)5), ((double)6)} ;
	size_t h_startset[1] = {0} ;
	size_t h_countset[1] = {LEN_H};
	if ( nc_put_vara(ncid, h_id, h_startset, h_countset, h_data) ) ERR;
    }

    /* Bug caused nc_close to fail with NC_EHDFERR (HDF Error) */
    if (nc_close(ncid)) ERR;

    /* Check file can be opened and read correctly */
    {   		
	int format;
	int ndims, nvars, ngatts, xdimid, nunlim;
	nc_type lat_type, h_type;
	int lat_rank, lat_natts, h_rank, h_natts;
	if (nc_open(FILENAME, NC_NOWRITE, &ncid)) ERR;
	if ( nc_inq_format(ncid, &format) ) ERR;
	if ( format != NC_FORMAT_NETCDF4_CLASSIC ) ERR;
	if ( nc_inq(ncid, &ndims, &nvars, &ngatts, &xdimid) ) ERR;
	if ( nc_inq_varid(ncid, NAME_LAT, &lat_id) ) ERR;
	if ( nc_inq_var(ncid, lat_id, NULL, &lat_type, &lat_rank, lat_dims, &lat_natts) ) ERR;
	if ( lat_type != NC_DOUBLE || lat_rank != RANK_LAT || lat_natts != 0 ) ERR;
	if ( nc_inq_varid(ncid, NAME_H, &h_id) ) ERR;
	if ( nc_inq_var(ncid, h_id, NULL, &h_type, &h_rank, h_dims, &h_natts) ) ERR;
	if ( h_type != NC_DOUBLE || h_rank != RANK_H || h_natts != 0 ) ERR;
	{
	    double lat_data[LEN_LAT];
	    size_t start[RANK_LAT] = {0} ;
	    size_t count[1] = {LEN_LAT};
	    if ( nc_get_vara(ncid, lat_id, start, count, lat_data) ) ERR;
	    if ( lat_data[0] != -45.0 || lat_data[1] != 45.0 ) ERR;
	}
	{
	    double h_data[LEN_H];
	    size_t start[RANK_H] = {0} ;
	    size_t count[1] = {LEN_H};
	    if ( nc_get_vara(ncid, h_id, start, count, h_data) ) ERR;
	    if ( h_data[0] != 5 || h_data[1] != 6 ) ERR;
	}
    }
    if (nc_close(ncid)) ERR;
      
   SUMMARIZE_ERR;
   FINAL_RESULTS;
}
