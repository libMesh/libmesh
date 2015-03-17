#!/bin/sh

# This script runs some parallel-netcdf (pnetcdf) I/O tests

set -e
echo
echo "Testing file created with pnetcdf is modifiable with netCDF..."
./tst_pnetcdf

set -x
# We assume a min of at least 2 processors is available
mpiexec -n 2 ./tst_parallel2
