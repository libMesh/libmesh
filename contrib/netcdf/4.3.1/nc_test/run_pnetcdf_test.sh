#!/bin/sh

# This script runs some parallel-netcdf (pnetcdf) I/O tests

set -e
echo
echo "Testing file created with pnetcdf is modifiable with netCDF..."
./tst_pnetcdf

set -x
mpiexec -n 4 ./tst_parallel2
