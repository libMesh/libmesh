#include <stdlib.h>
#include "netcdf.h"

int
main(void)
{
   int ncid, typeid, varid;
   float missing_value = -999.0;
   nc_vlen_t missing_val;

   if(nc_create("tst_vlen_data.nc", NC_CLOBBER | NC_NETCDF4, &ncid)) abort();

   /* Create a vlen type. */
   if (nc_def_vlen(ncid, "row_of_floats", NC_FLOAT, &typeid)) abort();;

   /* Declare a scalar variable of the vlen type */
   if (nc_def_var(ncid, "ragged_array", typeid, 0, NULL, &varid)) abort();;

   /* Create and write a variable attribute of the vlen type */
   missing_val.p = &missing_value;
   missing_val.len = 1;
   if (nc_put_att(ncid, varid, "_FillValue", typeid, 1, (void *) &missing_val)) abort();;

   if (nc_close(ncid)) abort();;
   exit(0);
}
