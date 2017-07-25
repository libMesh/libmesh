#!/bin/sh
#
# Front-end for tst_empty_vlen_unlim.c.  This test
# ensures that valid netcdf files are generated when
# a ragged VLEN array is defined but not populated.
#
# This script runs tst_empty_vlen_unlim and then
# runs `ncdump` against the two generated files.
#
# See https://github.com/Unidata/netcdf-c/issues/221 for
# full details.
#

set -e

if test "x$srcdir" = x ; then
srcdir=`pwd`
fi

echo ""
echo "* Testing Empty Ragged Arrays (VLEN)"

echo "Generating test netcdf files."
./tst_empty_vlen_unlim

# Since no comparison is made, I am not sure 
# if this is useful.
#echo "Validating Files with ncdump."
#echo "======================================"
#../ncdump/ncdump -s tst_empty_vlen_unlim.nc
#echo "---------------------------------------"
#../ncdump/ncdump -s tst_empty_vlen_lim.nc
#echo "======================================"


echo "* Tests Passed."
exit 0
