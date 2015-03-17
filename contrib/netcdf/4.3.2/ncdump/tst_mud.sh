#!/bin/sh
# This shell script tests ncdump and ncgen on netCDF-4 variables with multiple 
# unlimited dimensions.
# $Id $

set -e

if test "x$srcdir" = "x"; then
    srcdir=`dirname $0`; 
fi
# add hack for sunos
export srcdir;

echo ""
echo "*** Testing ncdump output for multiple unlimited dimensions"
echo "*** creating netcdf file tst_mud4.nc from ref_tst_mud4.cdl ..."
../ncgen/ncgen -k3 -b -o tst_mud4.nc $srcdir/ref_tst_mud4.cdl
echo "*** creating tst_mud4.cdl from tst_mud4.nc ..."
./ncdump tst_mud4.nc > tst_mud4.cdl
# echo "*** comparing tst_mud4.cdl with ref_tst_mud4.cdl..."
diff -b tst_mud4.cdl $srcdir/ref_tst_mud4.cdl
# echo "*** comparing annotation from ncdump -bc tst_mud4.nc with expected output..."
./ncdump -bc tst_mud4.nc > tst_mud4-bc.cdl
diff -b tst_mud4-bc.cdl $srcdir/ref_tst_mud4-bc.cdl
echo "*** All ncdump test output for multiple unlimited dimensions passed!"
exit 0
