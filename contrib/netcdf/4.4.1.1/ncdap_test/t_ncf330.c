/*! Test program for netcdf issue NCF-330
 *
 * This test was provided by Ellen Johnson at Mathworks and
 * illustrates an issue currently only seen on Windows.
 *
 * See https://bugtracking.unidata.ucar.edu/browse/NCF-330
 */

#include "config.h"
#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif
#ifdef HAVE_STDIO_H
#include <stdio.h>
#endif
#include <netcdf.h>

#ifdef _MSC_VER
#include <malloc.h>
#include <crtdbg.h>
#endif

static char* URL="http://data.nodc.noaa.gov/thredds/dodsC/testdata/pathfinderAgg/pathFinderV5.2_night.ncml";

int
main()
{
     int ncid;
     int ncstatus;
     int lat_id;
     nc_type lat_type;
     int lat_ndims;
     int lat_dimids[NC_MAX_VAR_DIMS];
     int lat_natts;
     int format_p;
     float lat_data[4320];

#ifdef _MSC_VER
	 /* See https://msdn.microsoft.com/en-us/library/5at7yxcs(v=vs.120).aspx
	    for more information re: these. */
	 /*int tmpFlag = _CrtSetDbgFlag(_CRTDBG_CHECK_ALWAYS_DF | _CRTDBG_LEAK_CHECK_DF);*/
#endif
     printf(" \n");
     printf("********************\n");
     printf("open URL %s\n",URL);
     printf(" \n");
     ncstatus = nc_open(URL, NC_NOWRITE, &ncid);
     if(ncstatus != NC_NOERR) {
         printf("Could not open: %s; server may be down; test ignored\n",URL);
         exit(0);
     }
     printf("status after open = %d\n", ncstatus);

     /* get the format */
     ncstatus = nc_inq_format(ncid, &format_p);
     printf("lat id = %d\n", format_p);

     /* get varid for latitude */
     ncstatus = nc_inq_varid(ncid,"lat",&lat_id);
     printf("status after inq lat id = %d\n", ncstatus);
     printf("lat id = %d\n", lat_id);

     ncstatus = nc_inq_var(ncid, lat_id, 0, &lat_type, &lat_ndims, lat_dimids, &lat_natts);
     printf("status after inq lat var = %d\n", ncstatus);
     printf("lat type = %d\n", lat_type);

     /* extract the first latitude value */
     ncstatus = nc_get_var_float(ncid, lat_id, &lat_data[0]);
     printf("status after get = %d\n", ncstatus);
     printf("my datum = %f\n", lat_data[0]);

     /* This code works okay in Linux and Mac
        Everything works okay up until here in Windows then an exception occurs
     */
     ncstatus = nc_close(ncid);
     printf("status after close = %d\n", ncstatus);
     printf("End of test.\n\n");

     return 0;
}
