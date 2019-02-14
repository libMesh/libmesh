#!/bin/sh

# This script runs some PnetCDF I/O tests

set -e
echo
echo "Testing file created with PnetCDF is modifiable with netCDF..."
./tst_pnetcdf

echo "Testing file created with PnetCDF works when adding variables..."
./tst_addvar tst_pnetcdf.nc

# We assume a min of at least 2 processors is available
mpiexec -n 2 ./tst_parallel2
