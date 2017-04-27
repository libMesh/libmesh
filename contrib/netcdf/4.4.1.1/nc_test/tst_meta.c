/*
  Copyright 20014, UCAR/Unidata
  See COPYRIGHT file for copying and redistribution conditions.

  This is part of netCDF.

  This program checks to see that netcdf_meta.h exists and is
  properly formatted.

*/

#include <stdio.h> /* printf() */
#include <config.h>
#include <stdlib.h>
#include <stdio.h>
#include <netcdf.h>
#include <netcdf_meta.h>

int main(int argc, char **argv) {

  /* If netcdf_meta.h is damaged, this file will
     just flat-out fail to compile, also resulting
     in an error.
  */

#ifndef NETCDF_META_H
  printf("Error! NETCDF_META_H not defined. Check netcdf_meta.h.\n");
  return -1;
#else
  printf("Success! NETCDF_META_H defined.\n");
  return 0;
#endif


}
