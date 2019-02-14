/* This is part of the netCDF package.
   Copyright 2016 University Corporation for Atmospheric Research/Unidata
   See COPYRIGHT file for conditions of use.

   Provided in support of https://github.com/Unidata/netcdf-c/issues/282
   Test provided by Greg Sjaardema

   Tests to see if the hashmap is being properly updated.

   */
#include <config.h>
#include <netcdf.h>

#define FILENAME "tst_rehash.nc"

int main()
{
  int  status;
  int  id;
  int  v1, v2, v3, v4;
  int  dimids[2];


  nc_create(FILENAME, NC_CLOBBER, &id);
  nc_redef(id);

  if ((status = nc_def_dim(id, "dim1", 10, &dimids[0])))
      return status;
  if ((status = nc_def_var(id, "dim1", NC_FLOAT, 1, dimids, &v1)))
      return status;
  if ((status = nc_def_var(id, "var1", NC_FLOAT, 1, dimids, &v2)))
      return status;

  nc_close(id);

  nc_open(FILENAME, NC_WRITE, &id);

  nc_redef(id);
  nc_rename_var(id, v1,"dim_new1");
  nc_rename_dim(id, dimids[0], "dim_new1");

  if ((status = nc_def_dim(id, "dim2", 20, &dimids[1])))
      return status;
  nc_def_var(id, "dim2", NC_FLOAT, 1, &dimids[1], &v3);
  nc_def_var(id, "var2", NC_FLOAT, 2, dimids,    &v4);

  nc_close(id);
  return 0;
}
